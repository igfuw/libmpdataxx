#include <cmath>
#include <boost/math/constants/constants.hpp>

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

#include "reversing_deform.hpp"

using namespace libmpdataxx;

struct ct_test_params_t
{
  // choosing non-dimensional sphere radius and evolution period
  static constexpr ::T
    R = 1,
    T = 5;
  // calculate advectors at nodal points and then average to staggered positions
  // or calculate at the staggered positions
  static constexpr bool node = false;
};
// to make refering to compile time test params easier
using tp = ct_test_params_t;

// gaussian hills initial condition
struct gauss_t
{
  T x0, y0;
  static constexpr T b = 5, hmax = 0.95;
  T xc(const T x, const T y) const { return tp::R * std::cos(y) * std::cos(x); }
  T yc(const T x, const T y) const { return tp::R * std::cos(y) * std::sin(x); }
  T zc(const T x, const T y) const { return tp::R * std::sin(y); }
  
  T operator()(const T x, const T y) const
  {
    return hmax * std::exp(-b * ( std::pow(xc(x, y) - xc(x0, y0), 2)
                                + std::pow(yc(x, y) - yc(x0, y0), 2)
                                + std::pow(zc(x, y) - zc(x0, y0), 2)));
  }
  BZ_DECLARE_FUNCTOR2(gauss_t);
};

// cosine bells initial condition
struct bells_t
{
  T x0, y0;
  static constexpr T r = tp::R / 2, hmax = 1.0, c = 0.9;
  T operator()(const T x, const T y) const
  {
    T ri = tp::R * std::acos(std::sin(y0) * std::sin(y) + std::cos(y0) * std::cos(y) * std::cos(x - x0));
    return ri < r ? c * hmax / 2 * (1 + std::cos(pi * ri / r)) : 0;
  }
  BZ_DECLARE_FUNCTOR2(bells_t);
};

// slotted cylinders initial condition
struct scyls_t
{
  T x0, y0;
  bool flip;
  static constexpr T r = tp::R / 2, c = 1.0;
  T operator()(const T x, const T y) const
  {
    T ri = tp::R * std::acos(std::sin(y0) * std::sin(y) + std::cos(y0) * std::cos(y) * std::cos(x - x0));
    bool cond1 = ri < r && std::abs(x - x0) >= r / (6 * tp::R);
    bool cond2 = flip ? (y - y0) < -5. / 12 * r / tp::R : (y - y0) > 5. / 12 * r / tp::R;
    bool cond3 = ri < r && (std::abs(x - x0) <  r / (6 * tp::R) && cond2);

    return cond1 || cond3 ? c : 0;
  }
  BZ_DECLARE_FUNCTOR2(scyls_t);
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
    enum { n_eqns = 4 };
    enum { opts = opts_arg };
    enum { sptl_intrp = tp::node ? solvers::exact : solvers::aver2};
  };
  

  const int nx = 2 * ny + 1;
  const T dx = 2 * pi / (nx - 1), dy = pi / ny;
  const T time = tp::T;
  const T dt = dx / (48 * 2 * pi);
  const int nt = time / dt;

  using slv_out_t = 
      output::hdf5_xdmf<
        reversing_deform<ct_params_t, ct_test_params_t>
    >;
  typename slv_out_t::rt_params_t p;

  p.n_iters = opts_iters; 
  p.grid_size = {nx, ny};
  p.dt = dt;
  p.di = dx;
  p.dj = dy;
  p.max_courant = max_cfl;

  p.outfreq = var_dt_arg ? tp::T / 2 : nt / 2; 
  p.outvars[0].name = "gh";
  p.outvars[1].name = "cb";
  p.outvars[2].name = "sc";
  p.outvars[3].name = "ccb";
  p.outdir = "out_" + dir_name;
  
  concurr::threads<
    slv_out_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::polar, bcond::polar
  > run(p); 

  blitz::firstIndex i;
  blitz::secondIndex j;

  // coordinates
  decltype(run.advectee()) X(run.advectee().extent()), Y(run.advectee().extent());
  X = i * dx;
  Y = (j + 0.5) * dy - pi / 2;

  std::array<T, 2> x0s = {5 * pi / 6, 7 * pi / 6};
  std::array<T, 2> y0s = {0, 0};

  // initial conditions
  const T b = 0.1;
  run.advectee(0) = 0;
  run.advectee(1) = b;
  run.advectee(2) = b;

  for (int m = 0; m < 2; ++m)
  {
    run.advectee(0) += gauss_t{x0s[m], y0s[m]}(X, Y);
    run.advectee(1) += bells_t{x0s[m], y0s[m]}(X, Y);
    run.advectee(2) += scyls_t{x0s[m], y0s[m], m == 0}(X, Y);
  }
  run.advectee(3) = -0.8 * blitz::pow(run.advectee(1), 2) + 0.9;
  run.g_factor() = tp::R * blitz::cos(Y) * dx * dy;
  
  // to avoid divergence check
  run.advector(0) = 0;
  run.advector(1) = 0;

  // integration
  run.advance(var_dt_arg ? time : nt);
}
