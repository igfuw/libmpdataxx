#pragma once

#include <cmath>

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include "moving_vort.hpp"

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
void test(const std::string &base_name, const int ny, const T max_cfl)
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
    enum { sptl_intrp = solvers::aver2 };
  };
  
  const int nx = 2 * ny + 1;
  const T dx = 2 * pi / (nx - 1), dy = pi / ny;
  const T time = 12;
  const T dt = dx / (48 * 2 * pi);
  const int nt = time / dt;

  using slv_out_t = 
      output::hdf5_xdmf<
        moving_vort<ct_params_t, ct_test_params_t>
    >;
  typename slv_out_t::rt_params_t p;

  p.n_iters = opts_iters; 
  p.grid_size = {nx, ny};
  p.dt = dt;
  p.di = dx;
  p.dj = dy;
  p.max_courant = max_cfl;

  p.outfreq = var_dt_arg ? 12.0 : nt;
  p.outvars[0].name = "psi";
  p.outdir = "out_" + dir_name;
  
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
}
