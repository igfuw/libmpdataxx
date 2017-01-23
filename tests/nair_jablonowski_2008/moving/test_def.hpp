/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <cmath>

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

#include "moving.hpp"
#include "../common/transforms.hpp"
#include "../common/convergence.hpp"

using namespace libmpdataxx;

struct ct_test_params_t
{
  // initial vortex position
  static constexpr T
    x0 = 3 * pi / 2,
    y0 = 0,
  // vortex velocity
    v0 = 2 * pi / 12,
  // solid-body rotation velocity
    u0 = 2 * pi / 12,
  // solid-body rotation angle
    a = pi / 2;
};

template <bool var_dt_arg, int opts_arg, int opts_iters>
stat_t<> test(const std::string &base_name, const int ny, const T max_cfl)
{
  auto dir_name = base_name + "_" + std::to_string(ny);
  std::cout << "Calculating: " << dir_name << std::endl;

  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { var_dt = var_dt_arg };
    enum { n_dims = 2 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };
  
  const int nx = 2 * ny + 1;
  const T dx = 2 * pi / (nx - 1), dy = pi / ny;
  const T time = 12;
  const T dt = dx / (48 * 2 * pi);
  const int nt = time / dt;

  using slv_out_t = 
      output::hdf5_xdmf<
        moving<ct_params_t, ct_test_params_t>
    >;
  typename slv_out_t::rt_params_t p;

  p.n_iters = opts_iters; 
  p.grid_size = {nx, ny};
  p.dt = dt;
  p.di = dx;
  p.dj = dy;
  p.max_courant = max_cfl;

  p.outfreq = var_dt_arg ? 6 : nt / 2; 
  p.outvars[0].name = "psi";
  p.outdir = dir_name;
  
  concurr::threads<
    slv_out_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::polar, bcond::polar
  > run(p); 
  
  // to make refering to compile time test params easier
  using tp = ct_test_params_t;
 
  // initialise coordinate transformations
  xpf_t xpf{.x0 = tp::x0, .y0 = tp::y0};
  ypf_t ypf{.x0 = tp::x0, .y0 = tp::y0};
  ixpf_t ixpf{.x0 = tp::x0, .y0 = tp::y0};
  iypf_t iypf{.x0 = tp::x0, .y0 = tp::y0};

  blitz::firstIndex i;
  blitz::secondIndex j;

  // coordinates
  decltype(run.advectee()) X(run.advectee().extent()), Y(run.advectee().extent());
  X = i * dx;
  Y = (j + 0.5) * dy - pi / 2;

  // helper arrays
  decltype(run.advectee()) r(run.advectee().extent()), omg(run.advectee().extent());

  // initial conditions
  r = 3 * cos(ypf(X, Y));
  omg = where(r != 0, tp::v0 * 3 * sqrt(2.) / (2 * r) * tanh(r) / pow2(cosh(r)), 0);

  run.advectee() = 1 - tanh(r / 5 * sin(xpf(X, Y)));
  run.g_factor() = blitz::cos(Y) * dx * dy;
  
  // to avoid divergence check
  run.advector(0) = 0;
  run.advector(1) = 0;

  // integration
  run.advance(var_dt_arg ? time : nt);

  // analytical solution
  decltype(run.advectee()) solution(run.advectee().extent());

  xpf.x0 = pi;
  xpf.y0 = pi / 2 - tp::a;
  ypf.x0 = pi;
  ypf.y0 = pi / 2 - tp::a;

  ixpf.x0 = pi;
  ixpf.y0 = pi / 2 - tp::a;
  iypf.x0 = pi;
  iypf.y0 = pi / 2 - tp::a;

  T x = xpf(tp::x0, tp::y0);
  T y = ypf(tp::x0, tp::y0);

  x += tp::u0 * run.time();

  T xc = ixpf(x, y);
  T yc = iypf(x, y);
  
  xpf.x0 = xc;
  xpf.y0 = yc;
  ypf.x0 = xc;
  ypf.y0 = yc;
  r = 3 * cos(ypf(X, Y));
  omg = where(r != 0, tp::v0 * 3 * sqrt(2.) / (2 * r) * tanh(r) / pow2(cosh(r)), 0);
  solution = 1 - tanh(r / 5 * sin(xpf(X, Y) - omg * run.time()));

  // calculation of stats
  stat_t<> ret;
  ret.min = min(run.advectee());
  ret.max = max(run.advectee());
  ret.L1 = sum(run.g_factor() * abs(run.advectee() - solution)) / sum(run.g_factor() * abs(solution));
  ret.L2 = sqrt(sum(run.g_factor() * pow2(run.advectee() - solution)) / sum(run.g_factor() * pow2(solution)));
  ret.Li = max(abs(run.advectee() - solution)) / max(abs(solution));
  return ret;
}
