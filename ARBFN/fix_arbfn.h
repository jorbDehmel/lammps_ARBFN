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
    Defines the `fix arbfn` class for extending LAMMPS. Based on FixQMMM
    from the QMMM package. Based on work funded by NSF grant
    2126451 at Colorado Mesa University.

    J Dehmel, J Schiffbauer, 2024/2025
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(arbfn,FixArbFn);
// clang-format on
#else

#ifndef FIX_ARBFN_HPP
#define FIX_ARBFN_HPP

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "interchange.h"

namespace LAMMPS_NS {
/**
 * @class FixArbFn
 * @brief A fix which communicates with an external controller
 * for arbitrary atomic force fixes. By default, this updates
 * every timestep. This is VERY slow.
 */
class FixArbFn : public Fix {
 public:
  FixArbFn(class LAMMPS *, int, char **);
  ~FixArbFn() override;

  void init() override;
  void post_force(int) override;
  int setmask() override;

 protected:
  uint controller_rank;
  double max_ms;
  MPI_Comm comm;
  uintmax_t every, counter;
  bool is_dipole = false;
  bool expect_response = true;
};
}    // namespace LAMMPS_NS

#endif    // FIX_ARBFN_HPP
#endif    // FIX_CLASS
