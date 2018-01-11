/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include "revolving_sphere_stats.hpp"
using namespace libmpdataxx;

template<int opts_arg, int opts_iters>
void test(const std::string& dir_name)
{
  enum {x, y, z};
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 3 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };

  int nt = 556;
  int nx = 59;

  using slv_t = solvers::mpdata<ct_params_t>;
  using slv_out_t = 
    stats<
      output::hdf5_xdmf<
        slv_t
      >
    >;
  typename slv_out_t::rt_params_t p;

  // pre instantation
  p.n_iters = opts_iters;
  p.grid_size = {nx, nx, nx};

  p.outfreq = nt;
  p.outvars[0].name = "psi";
  p.outdir = dir_name;

  // post instantation
  const typename ct_params_t::real_t
    dt = 0.018 * 2 * pi<double>(),
    L = 100,
    dx = L / (nx - 1),
    dy = dx,
    dz = dx,
    h = 4,
    r = 15,
    d = 25 / sqrt(3),
    x0 = 50 - d,
    y0 = 50 + d,
    z0 = 50 + d;
  
  p.di = dx;
  p.dj = dy;
  p.dk = dz;
  p.dt = dt;

  // instantation
  concurr::threads<
  slv_out_t,
  bcond::open, bcond::open,
  bcond::open, bcond::open,
  bcond::open, bcond::open
  > slv(p);

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;

  // sphere shape
  decltype(slv.advectee()) tmp(slv.advectee().extent());
  tmp.reindexSelf(slv.advectee().base());
  tmp =   blitz::pow(i * dx - x0, 2)
        + blitz::pow(j * dx - y0, 2)
        + blitz::pow(k * dx - z0, 2);
  slv.advectee() = where(tmp - pow(r, 2) <= 0, h, 0);

  const typename ct_params_t::real_t
    omega = 0.1,
    xc = 50,
    yc = 50,
    zc = 50;

  // constant angular velocity rotational field
  slv.advector(x) = omega / sqrt(3) * (-(j * dy - yc) + (k * dz - zc)) * dt / dx;
  slv.advector(y) = omega / sqrt(3) * ( (i * dx - xc) - (k * dz - zc)) * dt / dy;
  slv.advector(z) = omega / sqrt(3) * (-(i * dx - xc) + (j * dy - yc)) * dt / dz;

  // time stepping
  slv.advance(nt);
}

int main()
{
#if defined(USE_MPI)
  // we will instantiate many solvers, so we have to init mpi manually, 
  // because solvers will not know should they finalize mpi upon destruction
  MPI::Init_thread(MPI_THREAD_MULTIPLE);
#endif
  {
    enum { opts = 0 };
    enum { opts_iters = 1};
    test<opts, opts_iters>("upwind");
  }

  {
    enum { opts = opts::iga };
    enum { opts_iters = 2};
    test<opts, opts_iters>("iga");
  }

  {
    enum { opts = 0 };
    enum { opts_iters = 2};
    test<opts, opts_iters>("basic");
  }

  {
    enum { opts = opts::iga | opts::fct};
    enum { opts_iters = 2};
    test<opts, opts_iters>("iga_fct");
  }

  {
    enum { opts = opts::fct };
    enum { opts_iters = 2};
    test<opts, opts_iters>("fct");
  }
#if defined(USE_MPI)
  MPI::Finalize();
#endif
}
