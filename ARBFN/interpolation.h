/**
 * @file ARBFN/interpolation.h
 * @brief Provides interpolation functions for use in LAMMPS via
 * the ARBFN package. Specifically, these are for
 * `fix arbfn/ffield`.
 * @author J Dehmel, 2025. Written under MIT license.
 */

#pragma once

/**
 * @brief Linearly interpolates the 3-tuple of force deltas
 * between 2 points in space (by convention, we call the
 * positions x, but that doesn't need to be true).
 * @param _force_deltas Where to save the interpolated values
 * @param _x The distance from x0
 * @param _dx The the distance between x0 and x1
 * @param _x0_force_deltas The force deltas if you are on x0
 * @param _x1_force_deltas The force deltas if you are on x1
 */
inline void interpolate_line(double _force_deltas[3], const double &_x, const double &_dx,
                             const double _x0_force_deltas[3], const double _x1_force_deltas[3])
{
  const double x1_frac = (double) _x / (double) (_dx);
  const double x0_frac = 1.0 - x1_frac;
  _force_deltas[0] = _x0_force_deltas[0] * x0_frac + _x1_force_deltas[0] * x1_frac;
  _force_deltas[1] = _x0_force_deltas[1] * x0_frac + _x1_force_deltas[1] * x1_frac;
  _force_deltas[2] = _x0_force_deltas[2] * x0_frac + _x1_force_deltas[2] * x1_frac;
}

/**
 * @brief Bilinearly interpolate between 4 equidistant points
 * (by convention, we use x and y). _pos and _d_pos are
 * 3-tuples, but only the first 2 entries are used. It is
 * assumed that x1 > x0, y1 > y0
 * @param _force_deltas Where to save the interpolation
 * @param _pos The relative position to x0 y0
 * @param _d_pos The spacing between nodes (x, y, z)
 * @param _x0_y0_force_deltas The force deltas for x0y0
 * @param _x1_y0_force_deltas The force deltas for x1y0
 * @param _x0_y1_force_deltas The force deltas for x0y1
 * @param _x1_y1_force_deltas The force deltas for x1y1
 */
inline void interpolate_plane(double _force_deltas[3], const double _pos[3], const double _d_pos[3],
                              const double _x0_y0_force_deltas[3],
                              const double _x1_y0_force_deltas[3],
                              const double _x0_y1_force_deltas[3],
                              const double _x1_y1_force_deltas[3])
{
  double y0_force_deltas[3];
  double y1_force_deltas[3];

  interpolate_line(y0_force_deltas, _pos[0], _d_pos[0], _x0_y0_force_deltas, _x1_y0_force_deltas);
  interpolate_line(y1_force_deltas, _pos[0], _d_pos[0], _x0_y1_force_deltas, _x1_y1_force_deltas);

  interpolate_line(_force_deltas, _pos[1], _d_pos[1], y0_force_deltas, y1_force_deltas);
}

/**
 * @brief Trilinearly interpolate between 8 equidistant points
 * (by convention, x, y, and z). It is assumed that x1 > x0,
 * y1 > y0, z1 > z0.
 * @param _force_deltas The results of the interpolation
 * @param _pos The offset from x0y0z0
 * @param _d_pos The spacing between nodes
 * @param _x0_y0_z0_force_deltas The force deltas for x0 y0 z0
 * @param _x1_y0_z0_force_deltas The force deltas for x1 y0 z0
 * @param _x0_y1_z0_force_deltas The force deltas for x0 y1 z0
 * @param _x1_y1_z0_force_deltas The force deltas for x1 y1 z0
 * @param _x0_y0_z1_force_deltas The force deltas for x0 y0 z1
 * @param _x1_y0_z1_force_deltas The force deltas for x1 y0 z1
 * @param _x0_y1_z1_force_deltas The force deltas for x0 y1 z1
 * @param _x1_y1_z1_force_deltas The force deltas for x1 y1 z1
 */
inline void
interpolate_box(double _force_deltas[3], const double _pos[3], const double _d_pos[3],
                const double _x0_y0_z0_force_deltas[3], const double _x1_y0_z0_force_deltas[3],
                const double _x0_y1_z0_force_deltas[3], const double _x1_y1_z0_force_deltas[3],
                const double _x0_y0_z1_force_deltas[3], const double _x1_y0_z1_force_deltas[3],
                const double _x0_y1_z1_force_deltas[3], const double _x1_y1_z1_force_deltas[3])
{
  double z0_force_deltas[3];
  double z1_force_deltas[3];

  interpolate_plane(z0_force_deltas, _pos, _d_pos, _x0_y0_z0_force_deltas, _x1_y0_z0_force_deltas,
                    _x0_y1_z0_force_deltas, _x1_y1_z0_force_deltas);
  interpolate_plane(z1_force_deltas, _pos, _d_pos, _x0_y0_z1_force_deltas, _x1_y0_z1_force_deltas,
                    _x0_y1_z1_force_deltas, _x1_y1_z1_force_deltas);

  interpolate_line(_force_deltas, _pos[2], _d_pos[2], z0_force_deltas, z1_force_deltas);
}

/**
 * @brief Given some position, find the proper bin and
 * trilinearly interpolate with the 8 nearest points.
 * @param _force_deltas Where the results are stored
 * @param _pos The GLOBAL position (not relative!)
 * @param _minimal_pos The smallest edge of the simulation box
 * @param _nodes _nodes[x_index][y_index][z_index] is the
 * 3-tuple of force deltas at the given dimensional indices.
 * @param _position_deltas The "bin widths" between nodes
 * @param _num_nodes The count of nodes on each side
 */
inline void interpolate(double _force_deltas[3], const double _pos[3], const double _minimal_pos[3],
                        const double *const *const *const *const _nodes,
                        const double _position_deltas[3], const unsigned int _num_nodes[3])
{
  // Determine bin ids
  int x_bin = (int) ((_pos[0] - _minimal_pos[0]) / _position_deltas[0]);
  int y_bin = (int) ((_pos[1] - _minimal_pos[1]) / _position_deltas[1]);
  int z_bin = (int) ((_pos[2] - _minimal_pos[2]) / _position_deltas[2]);

  // Handle edge cases to avoid segfaults
  if (x_bin < 0) {
    x_bin = 0;
  } else if (x_bin + 1 >= _num_nodes[0]) {
    x_bin = _num_nodes[0] - 2;
  }
  if (y_bin < 0) {
    y_bin = 0;
  } else if (y_bin + 1 >= _num_nodes[1]) {
    y_bin = _num_nodes[1] - 2;
  }
  if (z_bin < 0) {
    z_bin = 0;
  } else if (z_bin + 1 >= _num_nodes[2]) {
    z_bin = _num_nodes[2] - 2;
  }

  // Localize to bins
  double local_position[3];
  local_position[0] = _pos[0] - (_minimal_pos[0] + x_bin * _position_deltas[0]);
  local_position[1] = _pos[1] - (_minimal_pos[1] + y_bin * _position_deltas[1]);
  local_position[2] = _pos[2] - (_minimal_pos[2] + z_bin * _position_deltas[2]);

  // Interpolate
  interpolate_box(_force_deltas, local_position, _position_deltas, _nodes[x_bin][y_bin][z_bin],
                  _nodes[x_bin + 1][y_bin][z_bin], _nodes[x_bin][y_bin + 1][z_bin],
                  _nodes[x_bin + 1][y_bin + 1][z_bin], _nodes[x_bin][y_bin][z_bin + 1],
                  _nodes[x_bin + 1][y_bin][z_bin + 1], _nodes[x_bin][y_bin + 1][z_bin + 1],
                  _nodes[x_bin + 1][y_bin + 1][z_bin + 1]);
}
