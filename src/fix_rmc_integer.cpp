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

#include "fix_rmc_integer.h"
#include "comm.h"
#include "modify.h"
#include "force.h"
#include "compute.h"
#include "pair.h"
#include "dihedral.h"
#include "bond.h"
#include "angle.h"
#include "improper.h"
#include "kspace.h"
#include "domain.h"
#include "update.h"
#include "random_park.h"
#include "error.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

FixRMCInteger::FixRMCInteger(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), c_pe(nullptr), random_equal(nullptr)
{ 
  c_pe = modify->get_compute_by_id("thermo_pe");
  
  //perform_step=0;
  perform_step=20000;  

  // First argument :- after how many MD steps to perform MC
  periodicity = utils::inumeric(FLERR, arg[3], false, lmp);
  
  // Second argument :- number of MC moves to perform per "turn"
  nmoves = utils::inumeric(FLERR, arg[4], false, lmp);

  // Third argument :- number of atoms in dopant molecule
  dopant_size= utils::inumeric(FLERR, arg[5], false, lmp);

  // Fourth argument :- number of atoms in semiconductor molecule
  semiconductor_size = utils::inumeric(FLERR, arg[6], false, lmp);

  // Fifth argument :- Number of dopant molecules
  n_dopant = utils::inumeric(FLERR, arg[7], false, lmp);

  // Sixth argument :- Number of semiconductor molecules
  n_semiconductor = utils::inumeric(FLERR, arg[8], false, lmp);

  // Seventh argument :- Temperature
  temperature = utils::numeric(FLERR, arg[9], false, lmp);

  // Eigth argument :- Reaction free energy
  delta_g = utils::numeric(FLERR, arg[10], false, lmp);

  // Ninth argument :- Type threshold (any atom type <= this is a semiconductor , >= is a dopant)
  // Critical variable to determine if a molecule is an semiconductor or a dopant
  type_threshold = utils::inumeric(FLERR, arg[11], false, lmp);


  num_semiconductor_charged=0;
  num_dopant_charged=0;
  num_semiconductor_neutral=n_semiconductor-num_semiconductor_charged;
  num_dopant_neutral=n_dopant-num_dopant_charged;
  
  beta = 1.0 / (force->boltz * temperature);
  
  acceptances=0;
  acceptance_rate=0; 
  rejections=0;
  nmcsteps=0;
  doping_efficiency=0; 
  neutral_to_charged_acceptances=0;
  neutral_to_charged_rejections=0;
 
  seed = 482794;
  random_equal = new RanPark(lmp, seed);

  size_limit = std::max(dopant_size, semiconductor_size);
 
  n_molecules=n_dopant+n_semiconductor;
  
  semiconductor_charged = new double [semiconductor_size];
  semiconductor_neutral = new double [semiconductor_size];
  dopant_charged = new double [dopant_size];
  dopant_neutral = new double [dopant_size];  

  TextFileReader semi_charge("semiconductor_charge.dat", "SCharges");
  TextFileReader semi_neutral("semiconductor_neutral.dat", "SNeutral");
  TextFileReader dope_charge("dopant_charge.dat", "DCharges");
  TextFileReader dope_neutral("dopant_neutral.dat", "DNeutral");
  semi_charge.next_dvector(semiconductor_charged, semiconductor_size);
  semi_neutral.next_dvector(semiconductor_neutral, semiconductor_size);
  dope_charge.next_dvector(dopant_charged, dopant_size);
  dope_neutral.next_dvector(dopant_neutral, dopant_size);

  
  //Read in the dihedrals that need to be altered
  TextFileReader dihedral_data("dihedral_list.dat", "dihedrals");
  
  //get number of dihedrals
  char *ndihedrals = dihedral_data.next_line(1);
  ValueTokenizer vt(ndihedrals, "\n");
  num_dihedrals = vt.next_int(); 
 
  
  //Get the dihedrals and the two types 
  //fmt::print(screen, "{} {}\n", "The number of dihedrals to process is ", num_dihedrals); 
  dihedral_list = new int*[num_dihedrals];
  for (int i=0;i<num_dihedrals;i++)
  {
     dihedral_list[i] = new int[7];
     char *dihedral_line = dihedral_data.next_line(6);
     ValueTokenizer vt(dihedral_line);
     for (int j=0;j<6;j++)
     {
        dihedral_list[i][j] = vt.next_int();
        fmt::print(screen, "{} ", dihedral_list[i][j]);
     }
     dihedral_list[i][6] = determine_molecule(dihedral_list[i][0]);
  }
  

  if (comm->me == 0)
  {
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen, "{}\n", "              RMC INITIALIZATION               ");
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen,"{} {}\n","RMC frequency: ", periodicity);
      fmt::print(screen,"{} {}\n","Number of RMC moves per turn: ",nmoves);
      fmt::print(screen,"{} {}\n","Size of dopant molecule: ",dopant_size);
      fmt::print(screen,"{} {}\n","Size of semiconductor molecule: ",semiconductor_size);
      fmt::print(screen,"{} {}\n","Number of dopant molecules: ",n_dopant);
      fmt::print(screen,"{} {}\n","Number of semiconductor molecules: ",n_semiconductor);
      fmt::print(screen,"{} {}\n","Temperature: ",temperature);
      fmt::print(screen,"{} {}\n","Reaction Energy: ", delta_g);
      fmt::print(screen,"{} {}\n","Type threshold: ",type_threshold);
      fmt::print(screen,"{} {}\n","betaDelta_G: ", beta*delta_g);
      fmt::print(screen, "{}\n", "###############################################");
  }
  MPI_Barrier(world);
}

int FixRMCInteger::determine_molecule(int global_id)
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

int FixRMCInteger::setmask()
{
  int mask = 0;
  mask |= FixConst::INITIAL_INTEGRATE;
  return mask; 
}

double FixRMCInteger::energy_full()
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

FixRMCInteger::Mol FixRMCInteger::initialize_molecule(int num_atoms_max)
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

int FixRMCInteger::determine_dopant_or_semiconductor(int mol_id)
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

int FixRMCInteger::determine_charged_or_neutral(struct Mol* molecule, int d_or_s)
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
          if (test_charge == semiconductor_charged[local_tag-1])
          {
             //fmt::print(screen, "{}\n", "Semiconductor is charged");
             charge_indicator[comm->me] = 1; // charged molecule, need to make neutral
          }
          else if (test_charge == semiconductor_neutral[local_tag-1])
          {
             //fmt::print(screen, "{}\n", "Semiconductor is neutral");
             charge_indicator[comm->me] = 0; // neutral molecule, need to make charged
          }
          else
          {
             error->all(FLERR, "No charge match for this semiconductor! Check the input file");
          }
       }
       else if (d_or_s == 1)
       {
          if (test_charge == dopant_charged[local_tag-1])
          {
             //fmt::print(screen, "{}\n", "Dopant is charged");
             charge_indicator[comm->me] = 1; // charged dopant, need to make neutral
          }
          else if (test_charge == dopant_neutral[local_tag-1])
          {
             //fmt::print(screen, "{}\n", "Dopant is neutral");
             charge_indicator[comm->me] = 0; // neutral dopant, need to make charged
          }      
          else
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

FixRMCInteger::Mol FixRMCInteger::get_molecule(int mol_id, int num_atoms_max) 
{

  Mol molecule = initialize_molecule(num_atoms_max);
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
  molecule.local_atoms = atom_counter;
  lmin_tag[comm->me] = mintag;
  lcount_across_ranks[comm->me] = atom_counter;
  MPI_Allreduce(lcount_across_ranks, gcount_across_ranks, comm->nprocs, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(lmin_tag, gmin_tag, comm->nprocs, MPI_INT, MPI_SUM, world);
  MPI_Barrier(world);
  //fmt::print("Determine overall min tag\n");
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
  //fmt::print(screen, "{}, {}\n","Min tag val is ", gmintagval);
  for (int i=0;i<atom_counter;i++)
  {
    molecule.local_tag[i] = molecule.global_tag[i] - gmintagval + 1;
  }
  //fmt::print("Everything fine till here\n");
  //fmt::print(screen, "{},{},{},{}\n", gmin_tag[0], gmin_tag[1], gmin_tag[2], gmin_tag[3]);
  return molecule;
}


void FixRMCInteger::modify_charge(struct Mol *molecule, double *charge_list)
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

void FixRMCInteger::restore_charge(struct Mol *molecule)
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


void FixRMCInteger::delete_molecule(struct Mol *molecule)
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

void FixRMCInteger::post_mortem()
{
   acceptance_rate = (double)acceptances/(double)nmcsteps;
   doping_efficiency = (double)num_dopant_charged/(double)n_dopant;
 
   if (comm->me == 0)
   {
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen, "{}\n", "              RMC OUTPUT SUMMARY               ");
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen,"{} {}\n", "Number of RMC moves: ", nmcsteps);
      fmt::print(screen,"{} {}\n", "Final charged dopant number: ", num_dopant_charged);
      fmt::print(screen,"{} {}\n", "Final charged semiconductor number: ", num_semiconductor_charged);
      fmt::print(screen,"{} {}\n", "Doping efficiency: ", doping_efficiency);  
      fmt::print(screen,"{} {}\n", "Number of neutral to charged acceptances: ", neutral_to_charged_acceptances);
      fmt::print(screen,"{} {}\n", "Final acceptances: ", acceptances);
      fmt::print(screen,"{} {}\n", "Number of neutral to charged rejections: ", neutral_to_charged_rejections);
      fmt::print(screen,"{} {}\n", "Final rejections: ", rejections);
      fmt::print(screen,"{} {}\n", "Acceptance rate: ",acceptance_rate);
      fmt::print(screen, "{}\n", "###############################################");
   }
}


double FixRMCInteger::change_dihedral_parameters(int molecule_id, int which_change)
{
   // if which_change == 0, this means changing type from undoped to doped
   // elseif which_change == 1, this means changing type from doped to undoped
   
   double pre_energy = energy_full();
   
   // Circle through the relevant dihedrals and see if any need to be modified.
   for (int d=0;d<num_dihedrals;d++)
   {
      if (dihedral_list[d][6] == molecule_id)
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
                   if (which_change == 0)
                   // switch to doped type
                   {
                      //fmt::print(screen, "{} {} {} {}\n", "Found the dihedral, switching from type", 
                      //           atom->dihedral_type[i][j], "to", dihedral_list[d][5]);
                      atom->dihedral_type[i][j] = dihedral_list[d][5];
                   }
                   else if (which_change == 1)
                   // switch to undoped type
                   {
                     // fmt::print(screen, "{} {} {} {}\n", "Found the dihedral, switching from type", 
                     //            atom->dihedral_type[i][j], "to", dihedral_list[d][4]);
                      atom->dihedral_type[i][j] = dihedral_list[d][4];
                   }
                }
             } 
           }
        }
      }
   }
   MPI_Barrier(world);

   // Update all the dihedral information in the neighbor lists
   // since that is what is used in the energy calculation
   neighbor->build_topology();

   // Now calculate the new energy
   double new_energy = energy_full();
   double energy_diff = new_energy - pre_energy;
   return energy_diff;
}


void FixRMCInteger::initial_integrate(int /*vflag*/)
{
  //if (perform_step == 0 || perform_step == update->ntimestep)
  if (perform_step == update->ntimestep)
  {
     double edihedral = change_dihedral_parameters(1,0);
     /*
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
     double dde = (double) num_dopant_charged/(double)n_dopant;
     if (comm->me == 0)
     {
       fmt::print(screen, "{} {}\n", "Dynamic Doping Efficiency: ",dde);
     }
     MPI_Barrier(world);  

     // Update when next to perform ReactiveMC
     perform_step = update->ntimestep + periodicity;
     */
  }
}

void FixRMCInteger::make_move()
{
    double reaction_energy = 0.0;
    double prefactor = 0.0;
    double transition_probability = 0.0;
    double edihedral = 0.0;   
 
    // Initialize random number generator
    std::random_device device;
    std::mt19937 rng(device());
    std::uniform_int_distribution<> atom_dist(1,n_molecules);
    
    // Calculate the energy before we do any mischief
    double starting_energy = energy_full(); 
    //fmt::print(screen, "{}\n", n_molecules);
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
    Mol semiconductor = get_molecule(rand_semi, size_limit);
 
    // Determine if chosen semiconductor is charged or not
    int is_semiconductor_charged = determine_charged_or_neutral(&semiconductor, 0);
    //fmt::print(screen, "{}, {}\n", "The semiconductor charge state is ", is_semiconductor_charged);
    
    // Find a dopant, whose charge state is the same as the semiconductor
    Mol dopant;
    
    indicator=-1;
    c_indicator=-1;

    while (indicator != 1 || c_indicator != is_semiconductor_charged){
       if (comm->me == 0)
       {
          rand_dope = atom_dist(rng);
       }
       MPI_Bcast(&rand_dope, 1, MPI_INT, 0, world);
       indicator = determine_dopant_or_semiconductor(rand_dope);
       if (indicator == 1)
       {
          dopant = get_molecule(rand_dope, size_limit);
          c_indicator = determine_charged_or_neutral(&dopant, 1);
       }
    }
 
    MPI_Barrier(world);
    MPI_Bcast(&rand_dope, 1, MPI_INT, 0, world);
    //fmt::print(screen, "{}, {}\n", "We have randomly chosen dopant ", rand_dope);
    
    // Retrieve information on the chosen dopant and populate it in this struct
    dopant = get_molecule(rand_dope, size_limit);
  
    // Determine charge of dopant (again) and just verify that our filtering worked
    int is_dopant_charged = determine_charged_or_neutral(&dopant, 1);
    if (is_semiconductor_charged != is_dopant_charged)
    {
       error->all(FLERR, "Charge states don't match!");
    }
    //else
    //{
    //   fmt::print("Charge states match!\n");
    //}

    // Modify charge - if semiconductor/dopant neutral make it charged and vice versa
    if (is_semiconductor_charged == 0)
    {

       //Change dihedral coefficients
       edihedral = change_dihedral_parameters(rand_semi, 0);
       
       
       modify_charge(&semiconductor, semiconductor_charged);
       modify_charge(&dopant, dopant_charged);

       reaction_energy = delta_g;
       prefactor = (double) (num_semiconductor_neutral*num_dopant_neutral)/(double)((num_semiconductor_charged+1)*(num_dopant_charged+1));
    }
    else
    {
       //Change dihedral coefficients
       edihedral = change_dihedral_parameters(rand_semi, 1);

       modify_charge(&semiconductor, semiconductor_neutral);
       modify_charge(&dopant, dopant_neutral);

       reaction_energy = -delta_g;
       prefactor = (num_semiconductor_charged*num_dopant_charged)/((num_semiconductor_neutral+1)*(num_dopant_neutral+1));
    }
    
    // Recalculate energy with new charges
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
    if (transition_probability > random_equal->uniform())
    {
       // Move accepted
       acceptances = acceptances+1;
       if (comm->me == 0)
       {
          fmt::print(screen, "{}\n", "MOVE ACCEPTED");
       }
       if (is_semiconductor_charged == 0)
       {
          num_semiconductor_charged = num_semiconductor_charged+1;
          num_dopant_charged = num_dopant_charged+1;
          num_semiconductor_neutral = num_semiconductor_neutral-1;
          num_dopant_neutral = num_dopant_neutral-1;
          neutral_to_charged_acceptances = neutral_to_charged_acceptances + 1;
       }
       else
       {
          num_semiconductor_charged = num_semiconductor_charged-1;
          num_dopant_charged = num_dopant_charged-1;
          num_semiconductor_neutral = num_semiconductor_neutral+1;
          num_dopant_neutral = num_dopant_neutral+1;
       }
    }
    else
    {
       // Move rejected
       rejections = rejections+1;
       if (is_semiconductor_charged == 0)
       {
          // Revert dihedral coefficients
          edihedral = change_dihedral_parameters(rand_semi, 1);
          neutral_to_charged_rejections = neutral_to_charged_rejections + 1;
       }
       else if (is_semiconductor_charged == 1)
       {
          // Revert dihedral coefficients
          edihedral = change_dihedral_parameters(rand_semi, 0);
       }
       if (comm->me == 0)
       { 
          fmt::print(screen, "{}\n", "MOVE REJECTED!!!!");
       }
       // Restore the charges to what they were
       restore_charge(&semiconductor);
       restore_charge(&dopant);
    }
    
    // The step is done, so we can delete the molecules
    delete_molecule(&semiconductor);
    delete_molecule(&dopant);
}

FixRMCInteger::~FixRMCInteger()
{
   post_mortem();

   // Some cleanup
   delete semiconductor_charged;
   delete semiconductor_neutral;
   delete dopant_charged;
   delete dopant_neutral;
   semiconductor_charged = nullptr;
   semiconductor_neutral = nullptr;
   dopant_charged = nullptr;
   dopant_neutral = nullptr;
}
