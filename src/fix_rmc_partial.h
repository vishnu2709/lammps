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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rmc_partial, FixRMCPartial)
// clang-format on
#else

#ifndef LMP_FIX_RMC_PARTIAL_H
#define LMP_FIX_RMC_PARTIAL_H

#include "fix.h"
#include "compute.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "text_file_reader.h"
#include <random>

namespace LAMMPS_NS {

class FixRMCPartial: public Fix {
   public:
   FixRMCPartial(class LAMMPS *, int, char **);
   ~FixRMCPartial();
   int setmask() override;
   int determine_molecule(int);
   int determine_dopant_or_semiconductor(int);
   double energy_full();
   double change_dihedral_parameters(int, int);
   void initial_integrate(int) override;
   void make_move();
   void post_mortem();

   struct Mol {
      double **pos;
      double *charge;
      int *type;
      int *global_tag;
      int *local_tag;
      int local_atoms;
      int charge_state;
   };

   int determine_charge_state(struct Mol*, int);
   void modify_charge(struct Mol*, double*);
   void restore_charge(struct Mol*);
   void delete_molecule(struct Mol*);
   struct Mol initialize_molecule(int);
   struct Mol get_molecule(int, int, int);

   private:
   int perform_step;
   int periodicity;
   int nmoves;
   int dopant_size;
   int semiconductor_size;
   int size_limit;
   int n_dopant;
   int n_semiconductor;
   int n_molecules;
   int type_threshold;
   int num_charge_states;
   int num_dihedrals;
   
   int acceptances;
   int rejections;
   int nmcsteps;
   int *num_dopant_charge;
   int *num_semiconductor_charge;

   double acceptance_rate;
   double temperature;
   double beta;
   double delta_g; 
   double deltaQ;



   int *dihedral_types;
   int **dihedral_list;
   double *delta_g_list;
   double *charges;
   double *dde;
   double *doping_efficiency;
   double **semiconductor_charges;
   double **dopant_charges;
   class Compute *c_pe;
};
};
#endif
#endif
