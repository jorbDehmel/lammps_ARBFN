#include "../ARBFN/interpolation.h"
#include <cassert>
#include <cmath>
#include <iostream>

void assert_approx_eq(const double &_l, const double &_r, const double &_eps = 0.01)
{
  if (fabs(_l - _r) > _eps) {
    std::cout << "Error! " << _l << " != " << _r << " to within " << _eps << '\n';
  }
  assert(fabs(_l - _r) <= _eps);
}

int main()
{
  double out[3];
  double pos[3] = {0.0, 0.0, 0.0};
  const double bin_deltas[3] = {10.0, 20.0, 100.0};
  double left_vals[3];
  double right_vals[3];
  double exp_val[3];

  // Points to interpolate
  double x0y0z0[3] = {0.0, 0.0, 100.0};
  double x1y0z0[3] = {100.0, 0.0, -100.0};
  double x0y1z0[3] = {50.0, 50.0, -50.0};
  double x1y1z0[3] = {-50.0, 100.0, 0.0};
  double x0y0z1[3] = {-100.0, 100.0, 0.0};
  double x1y0z1[3] = {-50.0, 25.0, 100.0};
  double x0y1z1[3] = {10.0, 1.0, 99.0};
  double x1y1z1[3] = {-12.0, 34.0, 56.0};

  pos[0] = 4.0;
  interpolate_line(out, pos[0], bin_deltas[0], x0y0z0, x1y0z0);
  assert_approx_eq(out[0], 40.0);
  assert_approx_eq(out[1], 0.0);
  assert_approx_eq(out[2], 20.0);

  pos[0] = 2.5;
  interpolate_line(out, pos[0], bin_deltas[0], x0y1z0, x1y1z0);
  assert_approx_eq(out[0], 25.0);
  assert_approx_eq(out[1], 62.5);
  assert_approx_eq(out[2], -37.5);

  pos[0] = 5.0;    // out of 10.0
  pos[1] = 5.0;    // out of 20.0
  pos[2] = 5.0;    // out of 100.0

  interpolate_line(left_vals, pos[0], bin_deltas[0], x0y0z0, x1y0z0);
  interpolate_line(right_vals, pos[0], bin_deltas[0], x0y1z0, x1y1z0);
  exp_val[0] = left_vals[0] + (right_vals[0] - left_vals[0]) / 4.0;
  exp_val[1] = left_vals[1] + (right_vals[1] - left_vals[1]) / 4.0;
  exp_val[2] = left_vals[2] + (right_vals[2] - left_vals[2]) / 4.0;

  interpolate_plane(out, pos, bin_deltas, x0y0z0, x1y0z0, x0y1z0, x1y1z0);
  assert_approx_eq(out[0], exp_val[0]);
  assert_approx_eq(out[1], exp_val[1]);
  assert_approx_eq(out[2], exp_val[2]);

  interpolate_plane(left_vals, pos, bin_deltas, x0y0z0, x1y0z0, x0y1z0, x1y1z0);
  interpolate_plane(right_vals, pos, bin_deltas, x0y0z1, x1y0z1, x0y1z1, x1y1z1);
  exp_val[0] = left_vals[0] + (right_vals[0] - left_vals[0]) / 20.0;
  exp_val[1] = left_vals[1] + (right_vals[1] - left_vals[1]) / 20.0;
  exp_val[2] = left_vals[2] + (right_vals[2] - left_vals[2]) / 20.0;

  interpolate_box(out, pos, bin_deltas, x0y0z0, x1y0z0, x0y1z0, x1y1z0, x0y0z1, x1y0z1, x0y1z1,
                  x1y1z1);
  assert_approx_eq(out[0], exp_val[0]);
  assert_approx_eq(out[1], exp_val[1]);
  assert_approx_eq(out[2], exp_val[2]);

  return 0;
}
