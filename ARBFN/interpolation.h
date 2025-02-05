/**
 * @file ARBFN/interpolation.h
 * @brief Provides interpolation functions for use in LAMMPS via
 * the ARBFN package. Specifically, these are for
 * `fix arbfn/ffield`.
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <sys/types.h>

inline void interpolate_line(double _force_deltas[3], const double &_x, const double &_dx,
                             const double _x0_force_deltas[3], const double _x1_force_deltas[3])
{
  const double x1_frac = (double) _x / (double) (_dx);
  const double x0_frac = 1.0 - x1_frac;
  _force_deltas[0] = _x0_force_deltas[0] * x0_frac + _x1_force_deltas[0] * x1_frac;
  _force_deltas[1] = _x0_force_deltas[1] * x0_frac + _x1_force_deltas[1] * x1_frac;
  _force_deltas[2] = _x0_force_deltas[2] * x0_frac + _x1_force_deltas[2] * x1_frac;
}

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

/*
_force_deltas[x][y][z] is a 3-tuple of force deltas
_position_deltas is a 3-tuple of (dx, dy, dz)
_num_bins is a 3-tuple of (max_x_bin, max_y_bin, max_z_bin)
*/
inline void interpolate(double _force_deltas[3], const double _pos[3], const double _minimal_pos[3],
                        const double *const *const *const *const _nodes,
                        const double _position_deltas[3], const uint _num_nodes[3])
{
  // Determine bin ids
  int x_bin = (int) ((_pos[0] - _minimal_pos[0]) / _position_deltas[0]);
  int y_bin = (int) ((_pos[1] - _minimal_pos[1]) / _position_deltas[1]);
  int z_bin = (int) ((_pos[2] - _minimal_pos[2]) / _position_deltas[2]);

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
