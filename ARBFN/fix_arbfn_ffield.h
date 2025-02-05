/* -*- c++ -*- ----------------------------------------------------------
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    https://www.lammps.org/, Sandia National Laboratories
    LAMMPS development team: developers@lammps.org

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.

    See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------
    Defines the `fix arbfn/ffield` class for extending LAMMPS. Based on
    FixQMMM from the QMMM package. Based on work funded by NSF grant
    2126451 at Colorado Mesa University.

    J Dehmel, J Schiffbauer, 2024/2025
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(arbfn/ffield,FixArbFnFField);
// clang-format on
#else

#ifndef FIX_ARBFN_FFIELD_HPP
#define FIX_ARBFN_FFIELD_HPP

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "interchange.h"

namespace LAMMPS_NS {
/**
 * @class FixArbFnFField
 * @brief Sets up a static force field of a finite number of
 * nodes upon initialization via MPI. Uses a modified ARBFN
 * protocol to communicate with the controller. Particles have
 * fixes applied to them by linear interpolation based on their
 * POSITIONS between nodes. This is independent of their
 * velocities or existing total force.
 */
class FixArbFnFField : public Fix {
 public:
  FixArbFnFField(class LAMMPS *, int, char **);
  ~FixArbFnFField() override;

  void init() override;
  void post_force(int) override;
  int setmask() override;

 protected:
  uint controller_rank;
  double max_ms;
  MPI_Comm comm;
  bool is_dipole = false;

  uint bin_counts[3];
  uint node_counts[3];

  double bin_deltas[3];
  double ****nodes = nullptr;
};
}    // namespace LAMMPS_NS

#endif    // FIX_ARBFN_FFIELD_HPP
#endif    // FIX_CLASS
