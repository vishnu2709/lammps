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
 *  ------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(rmc_integer, FixRMCInteger)
// clang-format on
#else

#ifndef LMP_FIX_RMC_INTEGER_H
#define LMP_FIX_RMC_INTEGER_H

#include "fix.h"
#include "compute.h"
#include "atom.h"
#include "text_file_reader.h"
#include <random>
namespace LAMMPS_NS {

class FixRMCInteger : public Fix {
   public:
   FixRMCInteger(class LAMMPS *, int, char **);
   ~FixRMCInteger() override;
   int setmask() override;
   void initial_integrate(int) override;
   struct Mol {
      double **pos;
      double *charge;
      int *type;
      int *global_tag;
      int *local_tag;
      int local_atoms;
   };

   
   int determine_charged_or_neutral(struct Mol*, int);
   struct Mol get_molecule(int, int);
   int determine_dopant_or_semiconductor(int);
   struct Mol initialize_molecule(int);
   void delete_molecule(struct Mol*);
   void modify_charge(struct Mol*, double*);
   void restore_charge(struct Mol*);
   void make_move();
   void post_mortem();
   int determine_molecule(int);
   double change_dihedral_parameters(int, int);


   double *semiconductor_charged;
   double *semiconductor_neutral;
   double *dopant_charged;
   double *dopant_neutral; 

   int num_dihedrals;
   int **dihedral_list;

   private:
   int perform_step, periodicity;
   int nmcsteps, nmoves;
   double energy_full();
   class Compute *c_pe;
   class RanPark *random_equal;
   int dopant_size, semiconductor_size; 
   int size_limit, type_threshold;
   int n_molecules, n_semiconductor, n_dopant;
   int num_semiconductor_charged, num_semiconductor_neutral;
   int num_dopant_charged, num_dopant_neutral;
   int acceptances, rejections;
   int neutral_to_charged_acceptances;
   int neutral_to_charged_rejections;
   double temperature, delta_g, beta;
   int seed;
   double acceptance_rate;
   double doping_efficiency;
};
};

#endif
#endif
