/**
 * @brief Defines library functions for use in the ARBFN library
 * @author J Dehmel, J Schiffbauer, 2024, MIT License
 */

#ifndef ARBFN_INTERCHANGE_H
#define ARBFN_INTERCHANGE_H

#include <list>
#include <mpi.h>

#define FIX_ARBFN_VERSION "0.2.1"

/**
 * @brief The color all ARBFN comms will be expected to have
 */
const static int ARBFN_MPI_COLOR = 56789;

/**
 * @struct AtomData
 * @brief Represents a single atom to be transferred
 */
struct AtomData {
  /// The x-position of the particle
  double x;

  /// The x-velocity of the particle
  double vx;

  /// The x-force of the particle
  double fx;

  /// The y-position of the particle
  double y;

  /// The y-velocity of the particle
  double vy;

  /// The y-force of the particle
  double fy;

  /// The z-position of the particle
  double z;

  /// The z-velocity of the particle
  double vz;

  /// The z-force of the particle
  double fz;

  /// Whether or not mu elements should be used
  bool is_dipole = false;

  /// x component of dipole moment orientation
  double mux;

  /// y component of dipole moment orientation
  double muy;

  /// z component of dipole moment orientation
  double muz;
};

/**
 * @struct FixData
 * @brief Represents a fix on an atom
 */
struct FixData {
  /// The delta to be added to fx
  double dfx;

  /// The delta to be added to fy
  double dfy;

  /// The delta to be added to fz
  double dfz;
};

/**
 * @struct FFieldNodeData
 * @brief Used for passing around FField fix nodes
 */
struct FFieldNodeData {
  /// The x index of the node
  unsigned int x_index;

  /// The y index of the node
  unsigned int ybin;

  /// The z index of the node
  unsigned int zbin;

  /// The delta to be added to fx
  double dfx;

  /// The delta to be added to fy
  double dfy;

  /// The delta to be added to fz
  double dfz;
};

/**
 * @brief Interchange, but for ffield fixes. This only happens once, upon simulation initialization.
 * @param _start A 3-tuple (x, y, z) of the lowest corner of the simulation box.
 * @param _bin_widths A 3-tuple for the x, y, and z spacing of the nodes.
 * @param _node_counts The number of nodes per side. A 3-tuple of the x, y, and z.
 * @param _controller_rank The rank of the controller within the provided communicator
 * @param _comm The MPI communicator to use
 * @returns A std::list of the data to be added
 */
std::list<FFieldNodeData> ffield_interchange(const double _start[3], const double _bin_widths[3],
                                             const unsigned int _bin_counts[3],
                                             const unsigned int &_controller_rank, MPI_Comm &_comm);

/**
 * @brief Send the given atom data, then receive the given fix data. This is blocking, but does not allow worker-side gridlocks.
 * @param _n The number of atoms/fixes in the arrays.
 * @param _from An array of atom data to send
 * @param _into An array of fix data that was received
 * @param _max_ms The max number of milliseconds to await each response
 * @param _controller_rank The rank of the controller within the provided communicator
 * @param _comm The MPI communicator to use
 * @returns true on success, false on failure
 */
bool interchange(const size_t &_n, const AtomData _from[], FixData _into[], const double &_max_ms,
                 const unsigned int &_controller_rank, MPI_Comm &_comm,
                 const bool &_expect_response = true);

/**
 * @brief Sends a registration packet to the controller.
 * @param _controller_rank The rank of the controller instance
 * @param _comm The communicator to use
 * @return True on success, false on error.
 */
bool send_registration(unsigned int &_controller_rank, MPI_Comm &_comm);

/**
 * @brief Sends a deregistration packet to the controller.
 * @param _controller_rank The MPI rank of the controller
 * @param _comm The communicator to use
 */
void send_deregistration(const int &_controller_rank, MPI_Comm &_comm);

#endif
