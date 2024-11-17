/* -*- c++ -*- ----------------------------------------------------------
 * LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 * https://www.lammps.org/, Sandia National Laboratories
 * LAMMPS development team: developers@lammps.org
 *          
 * Copyright (2003) Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.  This software is distributed under
 * the GNU General Public License.
 * 
 * See the README file in the top-level LAMMPS directory.
 * ------------------------------------------------------------------------- */

#include "fix_rmc_partial.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixRMCPartial::FixRMCPartial(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  perform_step = 0;

  // After how many MD steps should RMC be performed 
  periodicity = utils::inumeric(FLERR, arg[3], false, lmp);
  
  // Second argument :- number of MC moves to perform per "turn"
  nmoves = utils::inumeric(FLERR, arg[4], false, lmp);
  
  // Third argument :- number of atoms in dopant molecule
  dopant_size= utils::inumeric(FLERR, arg[5], false, lmp);
  
  // Fourth argument :- number of atoms in semiconductor molecule
  semiconductor_size = utils::inumeric(FLERR, arg[6], false, lmp);

  // Determine the max of the sizes, to use for memory allocation
  size_limit = std::max(dopant_size, semiconductor_size);
  
  // Fifth argument :- Number of dopant molecules
  n_dopant = utils::inumeric(FLERR, arg[7], false, lmp);
  
  // Sixth argument :- Number of semiconductor molecules
  n_semiconductor = utils::inumeric(FLERR, arg[8], false, lmp);
  
  // Defining the sum of them as a separate variable for convienience
  n_molecules = n_dopant+n_semiconductor;

  // Seventh argument :- Temperature
  temperature = utils::numeric(FLERR, arg[9], false, lmp);
  
  // Calculate beta
  beta = 1.0/(force->boltz * temperature);

  // Eigth argument :- Reaction free energy
  delta_g = utils::numeric(FLERR, arg[10], false, lmp);
  
  // Ninth argument :- Type threshold (any atom type <= this is a semiconductor , >= is a dopant)
  // Critical variable to determine if a molecule is an semiconductor or a dopant
  type_threshold = utils::inumeric(FLERR, arg[11], false, lmp);
  
  // Tenth argument (only for partial) :- How many charge states to study?
  // For rmc_integer, this is 2 by default (0 and 1)
  // But for rmc_partial, this can be specified by the user. 
  num_charge_states = utils::inumeric(FLERR, arg[12], false, lmp);


  // Separation between charge states
  deltaQ = 1.0/(num_charge_states-1.0);
  
  // List of charges
  charges = new double [num_charge_states];
  charges[0] = 0;
  for (int i=0;i<num_charge_states-1;i++)
  {
     charges[i+1] = charges[i] + deltaQ;
  }
  
  // Data structure to store semiconductor and dopant charges
  semiconductor_charges = new double*[num_charge_states];
  dopant_charges = new double*[num_charge_states];
  for (int i=0;i<num_charge_states;i++)
  {
     semiconductor_charges[i] = new double [semiconductor_size];
     dopant_charges[i] = new double [dopant_size];
  }

  // Initialize array to keep count of dopant/semiconductor in the different charge states
  num_dopant_charge = new int [num_charge_states];
  num_semiconductor_charge = new int[num_charge_states];
  num_dopant_charge[0] = n_dopant;
  num_semiconductor_charge[0] = n_semiconductor;
  for (int i=1;i<num_charge_states;i++)
  {
      num_dopant_charge[i] = 0;
      num_semiconductor_charge[i] = 0;
  }
  
  // Populate the neutral and integer charge state from files
  TextFileReader semi_integer("semiconductor_charge.dat", "SCharges");
  TextFileReader semi_neutral("semiconductor_neutral.dat", "SNeutral");
  TextFileReader dope_integer("dopant_charge.dat", "DCharges");
  TextFileReader dope_neutral("dopant_neutral.dat", "DNeutral");
  semi_neutral.next_dvector(semiconductor_charges[0], semiconductor_size);
  semi_integer.next_dvector(semiconductor_charges[num_charge_states-1], semiconductor_size);
  dope_neutral.next_dvector(dopant_charges[0], dopant_size);
  dope_integer.next_dvector(dopant_charges[num_charge_states-1], dopant_size);
  
  // Populate the intermediate charge values for each semiconductor atom using linear interpolation
  double temp_delta_semiconductor, temp_delta_dopant;
  
  for (int i=0;i<semiconductor_size;i++)
  {
     temp_delta_semiconductor = (semiconductor_charges[num_charge_states-1][i] - semiconductor_charges[0][i])/(num_charge_states-1);
     for (int j=0;j<num_charge_states-1;j++)
     {
        semiconductor_charges[j+1][i] = semiconductor_charges[j][i] + temp_delta_semiconductor;
     }
  }
  
  for (int i=0;i<dopant_size;i++)
  {
     temp_delta_dopant = (dopant_charges[num_charge_states-1][i] - dopant_charges[0][i])/(num_charge_states-1);
     for (int j=0;j<num_charge_states-1;j++)
     {
        dopant_charges[j+1][i] = dopant_charges[j][i] + temp_delta_dopant;
     }
  }
  
  // Create and populate the delta_g list for various charge states 
  delta_g_list = new double [num_charge_states];
  delta_g_list[0] = 0;
  double temp_delta_delta_g = delta_g/(num_charge_states-1.0);   
   
  for (int i=0;i<num_charge_states-1;i++)
  {
     delta_g_list[i+1] = delta_g_list[i] + temp_delta_delta_g;
  }

  // Read in the dihedrals that need to be altered
  TextFileReader dihedral_data("dihedral_list.dat", "dihedrals");

  // get number of dihedrals
  char *ndihedrals = dihedral_data.next_line(1);
  ValueTokenizer vt(ndihedrals, "\n");
  num_dihedrals = vt.next_int();

  // get dihedrals types
  dihedral_types = new int [num_charge_states];
  char *ndtypes = dihedral_data.next_line(num_charge_states);
  ValueTokenizer vt1(ndtypes);
  for (int i=0;i<num_charge_states;i++)
  {
     dihedral_types[i] = vt1.next_int();  
  }

  // Get the dihedral atom indices
  fmt::print(screen, "{} {}\n", "The number of dihedrals to process is ", num_dihedrals);  
  dihedral_list = new int*[num_dihedrals];
  
  for (int i=0;i<num_dihedrals;i++)
  {
     dihedral_list[i] = new int[5];
     char *dihedral_line = dihedral_data.next_line(4);
     ValueTokenizer vt(dihedral_line);
     for (int j=0;j<4;j++)
     {
        dihedral_list[i][j] = vt.next_int();
     }
     dihedral_list[i][4] = determine_molecule(dihedral_list[i][0]);
  }
   
   // Initialize the dynamic doping efficiency array
   dde = new double [num_charge_states];
   doping_efficiency = new double [num_charge_states];
   for (int i=0;i<num_charge_states;i++)
   {
      dde[i] = 0;
      doping_efficiency[i] = 0;
   }

   // Initialize acceptances/rejections
   acceptances = 0;
   rejections = 0;

   // get pointer for compute class, which will allow us to 
   // retrieve the potential energy
   c_pe = modify->get_compute_by_id("thermo_pe");

  if (comm->me == 0)
  {   
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen, "{}\n", "        RMC PARTIAL CHARGE INITIALIZATION      ");
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen,"{} {}\n","RMC frequency: ", periodicity);
      fmt::print(screen,"{} {}\n","Number of RMC moves per turn: ",nmoves);
      fmt::print(screen,"{} {}\n","Size of dopant molecule: ",dopant_size);
      fmt::print(screen,"{} {}\n","Size of semiconductor molecule: ",semiconductor_size);
      fmt::print(screen,"{} {}\n","Number of dopant molecules: ",n_dopant);
      fmt::print(screen,"{} {}\n","Number of semiconductor molecules: ",n_semiconductor);
      fmt::print(screen,"{} {}\n","Temperature: ",temperature);
      fmt::print(screen,"{} {}\n","Number of charge states: ", num_charge_states);
      fmt::print(screen,"{} ", "Charges:");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ",charges[i]);
      }
      fmt::print(screen,"\n{} ", "Reaction Energies: ");
      for (int i=0; i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ",delta_g_list[i]);
      }
      fmt::print(screen,"\n{} ", "betaDelta_G: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ", beta*delta_g_list[i]);
      }
      fmt::print(screen,"\n{} ", "Dihedral types: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ", dihedral_types[i]);
      }
      fmt::print(screen,"\n{} {}\n","Type threshold: ",type_threshold);
      fmt::print(screen, "{}\n", "###############################################");
  }
  MPI_Barrier(world);    
}

int FixRMCPartial::determine_molecule(int global_id)
{ 
  int mol_id = 0;
  int global_mol_id = 0;
  for (int i=0;i<atom->nlocal;i++)
  {  
     if (atom->tag[i] == global_id)
     { 
       mol_id = atom->molecule[i];
     }
  }
  MPI_Allreduce(&mol_id, &global_mol_id, 1, MPI_INT, MPI_SUM, world);
  return global_mol_id;
}

FixRMCPartial::Mol FixRMCPartial::initialize_molecule(int num_atoms_max)
{  
   Mol molecule;
   molecule.pos = new double*[num_atoms_max];
   for (int i=0;i<num_atoms_max;i++)
   { 
     molecule.pos[i] = new double[3];
   }
   molecule.charge = new double[num_atoms_max];
   molecule.type = new int[num_atoms_max];
   molecule.global_tag = new int[num_atoms_max];
   molecule.local_tag = new int[num_atoms_max];
   return molecule;
}

 
int FixRMCPartial::setmask()
{ 
  int mask = 0;
  mask |= FixConst::INITIAL_INTEGRATE;
  return mask;
}

void FixRMCPartial::initial_integrate(int /*vflag*/)
{ 
  if (perform_step == 0 || perform_step == update->ntimestep)
  {  
     for (int move=1;move<=nmoves;move++)
     {  
        if (comm->me == 0)
        {  
           fmt::print(screen, "{} {} {} {}\n", "Move ", move, " out of ", nmoves);
        }
        make_move();
        MPI_Barrier(world);
     }
     nmcsteps = nmcsteps + nmoves;
     
     // Calculate dynamic doping efficiency
     for (int i=0;i<num_charge_states;i++)
     {
         dde[i] = (double) num_dopant_charge[i]/(double)n_dopant;
     }
     

     if (comm->me == 0)
     {
        fmt::print(screen, "{}", "Dynamic Doping Efficiency: ");      
        for (int i=0;i<num_charge_states;i++)
        {
           fmt::print(screen, "{} ", dde[i]);
        }
        fmt::print(screen, "{}\n", " ");
     }
     MPI_Barrier(world);
     
     // Update when next to perform ReactiveMC
     perform_step = update->ntimestep + periodicity;
  }
}

double FixRMCPartial::energy_full()
{
  
  int eflag = 1;
  int vflag = 0;

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  if (modify->n_post_force_any) modify->post_force(vflag);
  
  double total_energy = c_pe->compute_scalar();
  update->eflag_global = update->ntimestep;

  return total_energy;
}

int FixRMCPartial::determine_dopant_or_semiconductor(int mol_id)
{
   int indicator[comm->nprocs];
   int global_indicator[comm->nprocs];
   int final_indicator;
   memset(indicator, 0.0, comm->nprocs*sizeof(int));
   memset(global_indicator, 0.0, comm->nprocs*sizeof(int));
   for (int i=0;i<atom->nlocal;i++)
   {
      indicator[comm->me] = -1;
      if (atom->molecule[i] == mol_id)
      {
        if (atom->type[i] > type_threshold)
        {
           indicator[comm->me] = 1; // dopant, above type_threshold
           break;
        }
        else
        {
           indicator[comm->me] = 0; // semiconductor, below type_threshold
           break;
        }
      }
   }
   MPI_Allreduce(indicator, global_indicator, comm->nprocs, MPI_INT, MPI_SUM, world);
   for (int procs=0;procs<comm->nprocs;procs++)
   {
      if (global_indicator[procs] != -1)
      {
         final_indicator = global_indicator[procs];
      }
   }
   return final_indicator;   
}

int FixRMCPartial::determine_charge_state(struct Mol* molecule, int d_or_s)
{
    int final_charge_indicator=-1;
    int charge_indicator[comm->nprocs];
    int global_charge_indicator[comm->nprocs];
    memset(charge_indicator, 0, comm->nprocs*sizeof(int));   
    memset(global_charge_indicator, 0, comm->nprocs*sizeof(int));
    if (molecule->local_atoms != 0)
    {
       double test_charge = molecule->charge[0];
       int local_tag = molecule->local_tag[0];
       //fmt::print(screen, "{}, {}, {}, {}\n", "The test charge is ", test_charge, " at local index ", local_tag); 
    
       if (d_or_s == 0)
       {
          charge_indicator[comm->me] = -1;
          for (int j=0;j<num_charge_states;j++)
          {
             if (test_charge == semiconductor_charges[j][local_tag-1])
             {
                charge_indicator[comm->me] = j;
             }
          }
          if (charge_indicator[comm->me] == -1)
          {
             error->all(FLERR, "No charge match for this semiconductor! Check the input file");
          }
       }
       else if (d_or_s == 1)
       {
          charge_indicator[comm->me] = -1;
          for (int j=0;j<num_charge_states;j++)
          {
            if (test_charge == dopant_charges[j][local_tag-1])
            {
               charge_indicator[comm->me] = j;
            }
          }
          if (charge_indicator[comm->me] == -1)
          {
             error->all(FLERR, "No charge match for this dopant! Check the input file");
          }
       }
    }
    else
    {
       charge_indicator[comm->me] = -1;
    }
    MPI_Barrier(world);
    MPI_Allreduce(charge_indicator, global_charge_indicator, comm->nprocs, MPI_INT, MPI_SUM, world);
    for (int procs=0;procs<comm->nprocs;procs++)
    {
       if (global_charge_indicator[procs] != -1)
       {
          final_charge_indicator = global_charge_indicator[procs];
       }
    }
    return final_charge_indicator;
}

double FixRMCPartial::change_dihedral_parameters(int molecule_id, int ending_state)
{
   // go to ending_state
   
   double pre_energy = energy_full();
   
   // Circle through the relevant dihedrals and see if any need to be modified.
   for (int d=0;d<num_dihedrals;d++)
   {
      if (dihedral_list[d][4] == molecule_id)
      { // this dihedral type needs to be modified
        // First we need to find the dihedral in the main data structure
        for (int i=0;i<atom->nlocal;i++) 
        {
           if (atom->tag[i] == dihedral_list[d][1])
           {
             for (int j=0;j<atom->num_dihedral[i];j++)
             {
                if (atom->dihedral_atom1[i][j] == dihedral_list[d][0] &&
                    atom->dihedral_atom2[i][j] == dihedral_list[d][1] &&
                    atom->dihedral_atom3[i][j] == dihedral_list[d][2] &&
                    atom->dihedral_atom4[i][j] == dihedral_list[d][3])
                {
                  // switch to new type
                  
                  //fmt::print(screen, "{} {} {} {}\n", "Found the dihedral, switching from type", 
                  //           atom->dihedral_type[i][j], "to", dihedral_types[ending_state]);
                  atom->dihedral_type[i][j] = dihedral_types[ending_state];
                }
             } 
           }
        }
      }
   }

   // Update all the dihedral information in the neighbor lists
   // since that is what is used in the energy calculation
   neighbor->build_topology();

   // Now calculate the new energy
   double new_energy = energy_full();
   double energy_diff = new_energy - pre_energy;
   return energy_diff;
}

FixRMCPartial::Mol FixRMCPartial::get_molecule(int mol_id, int num_atoms_max, int d_or_s) 
{
  // Initialize memory for the molecule
  Mol molecule = initialize_molecule(num_atoms_max);

  // Get all the atoms that belong to that molecule
  // This structure will look different across MPI ranks,
  // if atoms of a single molecule are split across ranks.
  int atom_counter=0; 
  int atom_counter_global=0;
  int lcount_across_ranks[comm->nprocs];
  int lmin_tag[comm->nprocs];
  int gmin_tag[comm->nprocs];
  int gcount_across_ranks[comm->nprocs];
  memset(lmin_tag, 0.0, comm->nprocs * sizeof(int));
  memset(lcount_across_ranks, 0.0, comm->nprocs * sizeof(int));
  int mintag = INT_MAX;
  for (int i=0;i< atom->nlocal;i++)
  {
    if (atom->molecule[i] == (double) mol_id)
    { 
      molecule.charge[atom_counter] = atom->q[i];
      molecule.type[atom_counter] = atom->type[i];
      molecule.global_tag[atom_counter] = atom->tag[i];
      molecule.local_tag[atom_counter] = 0.0;
    
      if (molecule.global_tag[atom_counter] < mintag)
      {
         mintag = molecule.global_tag[atom_counter];
      }
      for (int j=0;j<3;j++) 
      {
        molecule.pos[atom_counter][j] = atom->x[i][j];
      }
      atom_counter = atom_counter+1;
   }
  }
  // Add up all the atoms across the different MPI ranks
  molecule.local_atoms = atom_counter;
  lmin_tag[comm->me] = mintag;
  lcount_across_ranks[comm->me] = atom_counter;
  MPI_Allreduce(lcount_across_ranks, gcount_across_ranks, comm->nprocs, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(lmin_tag, gmin_tag, comm->nprocs, MPI_INT, MPI_SUM, world);
  MPI_Barrier(world);

  // Determine overall minimum global index
  int gmintagval = INT_MAX;
  for (int i=0;i<comm->nprocs;i++)
  {
    atom_counter_global = atom_counter_global + gcount_across_ranks[i];
    if (gmin_tag[i] != INT_MAX)
    {

      if (gmin_tag[i] < gmintagval)
      {
        gmintagval = gmin_tag[i];
      }
    }
  }
  // Subtract off the minimum global index so we get a 
  // new "local" index. This is not local to processor, but local to a molecule
  // This will be very helpful for charge manipulation
  for (int i=0;i<atom_counter;i++)
  {
    molecule.local_tag[i] = molecule.global_tag[i] - gmintagval + 1;
  }

  // Now, arguably the most important step
  // Determine the charge state of the molecule

  int charge_state = determine_charge_state(&molecule, d_or_s);
  molecule.charge_state = charge_state;
  return molecule;
}

void FixRMCPartial::modify_charge(struct Mol *molecule, double *charge_list)
{
  
  for (int i=0;i<molecule->local_atoms;i++)
  {
     for (int j=0;j<atom->nlocal;j++)
     {
        if (atom->tag[j] == molecule->global_tag[i])
        {
           atom->q[j] = charge_list[molecule->local_tag[i]-1];
        }
     }
  }
}

void FixRMCPartial::restore_charge(struct Mol *molecule)
{
  for (int i=0;i<molecule->local_atoms;i++)
  {
     for (int j=0;j<atom->nlocal;j++)
     {
        if (atom->tag[j] == molecule->global_tag[i])
        {
           atom->q[j] = molecule->charge[i];
        }
     }
  }

}

void FixRMCPartial::make_move()
{
    double reaction_energy = 0.0;
    double prefactor = 0.0;
    double transition_probability = 0.0;
    double edihedral = 0.0;
 
    // Initialize random number generator for atoms
    std::random_device device;
    std::mt19937 rng(device());
    std::uniform_int_distribution<> atom_dist(1,n_molecules);
    std::uniform_int_distribution<> type_dist(0,num_charge_states-1);
    std::uniform_real_distribution<> acc_dist(0, 1);

    
    // Calculate the energy before we do any mischief
    double starting_energy = energy_full(); 
    
    if (comm->me == 0)
    {
       fmt::print(screen, "{}, {}\n", "The starting energy is ", starting_energy);
    }
    int rand_semi, rand_dope;
    int indicator=-1;
    int c_indicator=-1;    
    
   
    // Find a semiconductor
    while (indicator != 0) {
       if (comm->me == 0)
       { 
         rand_semi = atom_dist(rng);
       }
       MPI_Bcast(&rand_semi, 1, MPI_INT, 0, world);
       indicator = determine_dopant_or_semiconductor(rand_semi);
    }
    
    MPI_Barrier(world);
    //fmt::print(screen, "{}, {}\n", "We have randomly chosen semiconductor ", rand_semi);

    // Retrieve information on the chosen molecule and populate it in this struct
    Mol semiconductor = get_molecule(rand_semi, size_limit, 0);


   // fmt::print(screen, "{}, {}\n", "The semiconductor charge state is ", charges[semiconductor.charge_state]);
    
    // Find a dopant, whose charge state is the same as the semiconductor
    Mol dopant;
    
    indicator=-1;
    c_indicator=-1;

    while (indicator != 1 || c_indicator != semiconductor.charge_state){
       if (comm->me == 0)
       {
          rand_dope = atom_dist(rng);
       }
       MPI_Bcast(&rand_dope, 1, MPI_INT, 0, world);
       indicator = determine_dopant_or_semiconductor(rand_dope);
       if (indicator == 1)
       {
          dopant = get_molecule(rand_dope, size_limit, 1);
          c_indicator = dopant.charge_state;
       }
    }

    //MPI_Bcast(&rand_dope, 1, MPI_INT, 0, world);
    MPI_Barrier(world);
    
    //fmt::print(screen, "{}, {}\n", "We have randomly chosen dopant ", rand_dope);
    
    // Retrieve information on the chosen dopant and populate it in this struct
    dopant = get_molecule(rand_dope, size_limit, 1);
  
    // Verify charge states are the same
    if (dopant.charge_state != semiconductor.charge_state)
    {
       error->all(FLERR, "Charge states don't match!");
    }

    // Identify a destination charge state, picked randomly
    // of course it has to be different from the starting state
    int destination_charge_state = dopant.charge_state;

    if (comm->me == 0)
    {
      while (destination_charge_state == dopant.charge_state)
      {
         destination_charge_state = type_dist(rng);
      }
    }
    MPI_Barrier(world);
    MPI_Bcast(&destination_charge_state, 1, MPI_INT, 0, world);
    if (comm->me == 0)
    {
      fmt::print(screen, "{} {} {} {}\n", "Going from charge state", charges[semiconductor.charge_state], "to ", charges[destination_charge_state]);
    }

    // Change the dihedral parameters to destination type 
    // and capture the dihedral energy
    edihedral = change_dihedral_parameters(rand_semi, destination_charge_state);

    // Modify charge to new type
    modify_charge(&semiconductor, semiconductor_charges[destination_charge_state]);
    modify_charge(&dopant, dopant_charges[destination_charge_state]);
    reaction_energy = delta_g_list[destination_charge_state] - delta_g_list[semiconductor.charge_state];
    prefactor = (double)(num_semiconductor_charge[semiconductor.charge_state]*num_dopant_charge[dopant.charge_state])/
      (double)(num_semiconductor_charge[destination_charge_state]+1)*(num_dopant_charge[destination_charge_state]+1);

    // Recalculate energy after these changes
    // subtract the dihedral energy, that has to be put back in with more thought later
    double new_energy = energy_full();
    double energy_diff = new_energy - starting_energy + reaction_energy - edihedral;

    if (comm->me == 0)
    {
       fmt::print(screen, "{} {} {} {}\n", "The new energy is ", new_energy, " and the difference is ", energy_diff);
       fmt::print(screen, "{} {} {}\n", "The dihedral energy difference is ", edihedral, " which was taken out of the acceptance criteria");
    }

    // Calculate acceptance probability
    transition_probability = prefactor*exp(-beta*energy_diff);
    //fmt::print(screen, "{}, {}\n", "The transition probability is ", transition_probability);
    
    // Determine whether to accept or reject
    if (transition_probability > acc_dist(rng))
    {
       // Move accepted
       acceptances = acceptances+1;
       if (comm->me == 0)
       {
          fmt::print(screen, "{}\n", "MOVE ACCEPTED");
       }
       MPI_Barrier(world);
       num_semiconductor_charge[destination_charge_state] = num_semiconductor_charge[destination_charge_state]+1;
       num_dopant_charge[destination_charge_state] = num_dopant_charge[destination_charge_state]+1;
       num_semiconductor_charge[semiconductor.charge_state] = num_semiconductor_charge[semiconductor.charge_state]-1;
       num_dopant_charge[dopant.charge_state] = num_dopant_charge[dopant.charge_state]-1;
    }
    else
    {
       // Move rejected
       rejections = rejections+1;

       // Revert dihedral coefficients
       edihedral = change_dihedral_parameters(rand_semi, semiconductor.charge_state);

       if (comm->me == 0)
       { 
          fmt::print(screen, "{}\n", "MOVE REJECTED!!!!");
       }
       MPI_Barrier(world);
       // Restore the charges to what they were
       restore_charge(&semiconductor);
       restore_charge(&dopant);
    }
    
    // The step is done, so we can delete the molecules
    delete_molecule(&semiconductor);
    delete_molecule(&dopant);
}

void FixRMCPartial::delete_molecule(struct Mol *molecule)
{
    delete molecule->type;
    delete molecule->charge;
    delete molecule->global_tag;
    delete molecule->local_tag;
    for (int i=0;i<molecule->local_atoms;i++)
    {
       delete molecule->pos[i];
    } 
    delete molecule->pos;
    molecule = nullptr;
}

void FixRMCPartial::post_mortem()
{
   acceptance_rate = (double)acceptances/(double)nmcsteps;
   double total_charged_dopants=0;
   double total_charged_semiconductors=0;
   for (int i=0;i<num_charge_states;i++)
   {
      doping_efficiency[i] = (double)num_dopant_charge[i]/(double)n_dopant;
      total_charged_dopants = total_charged_dopants+num_dopant_charge[i];
      total_charged_semiconductors = total_charged_semiconductors+num_semiconductor_charge[i];
   
   }
   double total_doping_efficiency = (double)total_charged_dopants/(double)n_dopant;
 
   if (comm->me == 0)
   {
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen, "{}\n", "              RMC OUTPUT SUMMARY               ");
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen,"{} {}\n", "Number of RMC moves: ", nmcsteps);
      fmt::print(screen, "{}", "Final charged dopant number: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen, "{} ", num_dopant_charge[i]);
      }
      fmt::print(screen,"\n{} {}\n", "Final charged dopant number: ", total_charged_dopants);
      fmt::print(screen,"{} {}\n", "Final charged semiconductor number: ", total_charged_semiconductors);
      fmt::print(screen,"{}", "Doping efficiency: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ", doping_efficiency[i]);
      }
      fmt::print(screen,"\n{} {}\n", "Final acceptances: ", acceptances);
      fmt::print(screen,"{} {}\n", "Final rejections: ", rejections);
      fmt::print(screen,"{} {}\n", "Acceptance rate: ",acceptance_rate);
      fmt::print(screen, "{}\n", "###############################################");
   }
}


FixRMCPartial::~FixRMCPartial()
{
   post_mortem();
}
