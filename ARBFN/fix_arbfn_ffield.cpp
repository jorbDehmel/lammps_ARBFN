#include "fix_arbfn_ffield.h"
#include "interchange.h"
#include "interpolation.h"
#include "utils.h"
#include <domain.h>
#include <mpi.h>
#include <string>

LAMMPS_NS::FixArbFnFField::FixArbFnFField(class LAMMPS *_lmp, int _c, char **_v) : Fix(_lmp, _c, _v)
{
  // Split comm
  PMPI_Comm_split(MPI_COMM_WORLD, ARBFN_MPI_COLOR, 0, &comm);

  // Handle keywords here
  max_ms = 0.0;

  if (_c < 6) {
    error->universe_one(FLERR, "Malformed `fix arbfn/ffield': Missing x/y/z bin counts.");
  }

  double bin_counts[3];
  bin_counts[0] = utils::numeric(FLERR, _v[3], false, _lmp);
  bin_counts[1] = utils::numeric(FLERR, _v[4], false, _lmp);
  bin_counts[2] = utils::numeric(FLERR, _v[5], false, _lmp);

  node_counts[0] = bin_counts[0] + 1;
  node_counts[1] = bin_counts[1] + 1;
  node_counts[2] = bin_counts[2] + 1;

  bin_deltas[0] = (double) (lmp->domain->boxhi[0] - lmp->domain->boxlo[0]) / (double) bin_counts[0];
  bin_deltas[1] = (double) (lmp->domain->boxhi[1] - lmp->domain->boxlo[1]) / (double) bin_counts[1];
  bin_deltas[2] = (double) (lmp->domain->boxhi[2] - lmp->domain->boxlo[2]) / (double) bin_counts[2];

  for (int i = 6; i < _c; ++i) {
    const char *const arg = _v[i];

    if (strcmp(arg, "every") == 0) {
      if (i + 1 >= _c) {
        error->universe_one(FLERR, "Malformed `fix arbfn': Missing argument for `every'.");
      }
      every = utils::numeric(FLERR, _v[i + 1], false, _lmp);
      ++i;
    } else if (strcmp(arg, "dipole") == 0) {
      is_dipole = !is_dipole;
    }

    else {
      error->universe_one(
          FLERR, "Malformed `fix arbfn/ffield': Unknown keyword `" + std::string(arg) + "'.");
    }
  }

  // Dynalloc
  nodes = new double ***[node_counts[0]];
  for (uint x_bin = 0; x_bin < node_counts[0]; ++x_bin) {
    nodes[x_bin] = new double **[node_counts[1]];
    for (uint y_bin = 0; y_bin < node_counts[1]; ++y_bin) {
      nodes[x_bin][y_bin] = new double *[node_counts[2]];
      for (uint z_bin = 0; z_bin < node_counts[2]; ++z_bin) {
        nodes[x_bin][y_bin][z_bin] = new double[3];
        nodes[x_bin][y_bin][z_bin][0] = 0.0;
        nodes[x_bin][y_bin][z_bin][1] = 0.0;
        nodes[x_bin][y_bin][z_bin][2] = 0.0;
      }
    }
  }
}

LAMMPS_NS::FixArbFnFField::~FixArbFnFField()
{
  send_deregistration(controller_rank, comm);
  MPI_Comm_free(&comm);

  if (nodes == nullptr) { return; }

  // Free bins here
  for (uint x_bin = 0; x_bin < node_counts[0]; ++x_bin) {
    for (uint y_bin = 0; y_bin < node_counts[1]; ++y_bin) {
      for (uint z_bin = 0; z_bin < node_counts[2]; ++z_bin) { delete[] nodes[x_bin][y_bin][z_bin]; }
      delete[] nodes[x_bin][y_bin];
    }
    delete[] nodes[x_bin];
  }
  delete[] nodes;
}

void LAMMPS_NS::FixArbFnFField::init()
{
  // Allocate bins here
  // double ****_all_force_deltas;
  //        |||^3-tuple for force deltas
  //        ||^z positions
  //        |^y positions
  //        ^x positions

  bool res = send_registration(controller_rank, comm);
  if (!res) {
    error->universe_one(
        FLERR, "`fix arbfn/ffield' failed to register with controller: Ensure it is running.");
  }

  // Populate bins from controller here
  // This is the first one, so we don't send any atomic data
  const auto points =
      ffield_interchange(lmp->domain->boxlo, bin_deltas, node_counts, controller_rank, comm, every);

  for (const auto &p : points) {
    if (p.x_index >= node_counts[0]) {
      error->universe_one(FLERR, "`fix arbfn/ffield' controller sent invalid x bin.");
    } else if (p.ybin >= node_counts[1]) {
      error->universe_one(FLERR, "`fix arbfn/ffield' controller sent invalid y bin.");
    } else if (p.zbin >= node_counts[2]) {
      error->universe_one(FLERR, "`fix arbfn/ffield' controller sent invalid z bin.");
    }
    nodes[p.x_index][p.ybin][p.zbin][0] += p.dfx;
    nodes[p.x_index][p.ybin][p.zbin][1] += p.dfy;
    nodes[p.x_index][p.ybin][p.zbin][2] += p.dfz;
  }
}

void LAMMPS_NS::FixArbFnFField::post_force(int)
{
  // Special refresh case
  if (every && ++counter >= every) {
    counter = 0;

    const double *const *const x = atom->x;
    const double *const *const v = atom->v;
    const double *const *const mu = atom->mu;
    const int *const mask = atom->mask;

    double *const *const f = atom->f;

    // Move from LAMMPS atom format to AtomData struct
    std::vector<AtomData> to_send;
    size_t n = 0;
    for (size_t i = 0; i < atom->nlocal; ++i) {
      if (mask[i] & groupbit) {
        AtomData to_add;
        to_add.x = x[i][0];
        to_add.y = x[i][1];
        to_add.z = x[i][2];
        to_add.vx = v[i][0];
        to_add.vy = v[i][1];
        to_add.vz = v[i][2];
        to_add.fx = f[i][0];
        to_add.fy = f[i][1];
        to_add.fz = f[i][2];

        to_add.is_dipole = is_dipole;
        if (to_add.is_dipole) {
          to_add.mux = mu[i][0];
          to_add.muy = mu[i][1];
          to_add.muz = mu[i][2];
        }

        to_send.push_back(to_add);
        ++n;
      }
    }

    const auto points =
        ffield_interchange(lmp->domain->boxlo, bin_deltas, node_counts, controller_rank, comm,
                           every, to_send.size(), to_send.data());

    for (const auto &p : points) {
      if (p.x_index >= node_counts[0]) {
        error->universe_one(FLERR, "`fix arbfn/ffield' controller sent invalid x bin.");
      } else if (p.ybin >= node_counts[1]) {
        error->universe_one(FLERR, "`fix arbfn/ffield' controller sent invalid y bin.");
      } else if (p.zbin >= node_counts[2]) {
        error->universe_one(FLERR, "`fix arbfn/ffield' controller sent invalid z bin.");
      }
      nodes[p.x_index][p.ybin][p.zbin][0] += p.dfx;
      nodes[p.x_index][p.ybin][p.zbin][1] += p.dfy;
      nodes[p.x_index][p.ybin][p.zbin][2] += p.dfz;
    }
  }

  const double *const *const x = atom->x;
  const int *const mask = atom->mask;
  double *const *const f = atom->f;
  double force_deltas[3];

  for (size_t i = 0; i < atom->nlocal; ++i) {
    if (mask[i] & groupbit) {
      interpolate(force_deltas, x[i], lmp->domain->boxlo, nodes, bin_deltas, node_counts);
      f[i][0] += force_deltas[0];
      f[i][1] += force_deltas[1];
      f[i][2] += force_deltas[2];
    }
  }
}

int LAMMPS_NS::FixArbFnFField::setmask()
{
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}
