/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

using T = double;
using namespace libmpdataxx;

template <int opts_arg, int opts_iters>
void test(const std::string filename)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 3 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };

  int nx = 257;
  int ny = 65;
  int nz = 65;

  T pi = boost::math::constants::pi<T>();
  T dx = 4.0 / (nx - 1);
  T dy = 1.0 / (ny - 1);
  T dz = 1.0 / (nz - 1);
  T dt = 0.8 * dx;

  T time = 4.0 + 3.0;
  int nt = time / dt;

  using slv_out_t = 
      output::hdf5_xdmf<
        solvers::mpdata<ct_params_t>
    >;
  typename slv_out_t::rt_params_t p;
  p.dt = dt;

  // pre instantiation
  p.n_iters = opts_iters; 
  p.grid_size = {nx, ny, nz};

  p.outfreq = nt / 8; 
  p.outvars[0].name = "psi";
  p.outdir = "out_3d";

  // instantiation
  concurr::threads<
    slv_out_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p); 

  T rr = 0.5;
  T x0 = 0.5;
  T y0 = 0.5;
  T z0 = 0.5;

  decltype(run.advectee()) r_local(run.advectee().extent());
  r_local.reindexSelf(run.advectee().lbound());

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;

  r_local = sqrt(blitz::pow(i * dx - x0, 2) + blitz::pow(j * dy - y0, 2) + blitz::pow(k * dz - z0, 2));
  run.advectee() = where(
    r_local <= rr,
    1 + 0.25 * pow(1 + cos(pi * r_local / rr), 2),
    1.
  );

  decltype(run.advectee()) solution(run.advectee_global().shape());
  decltype(run.advectee()) r_global(run.advectee_global().shape());

  T disp = std::fmod(time, 4.0);
  
  r_global = sqrt(blitz::pow(i * dx - (x0 + disp), 2) + blitz::pow(j * dy - y0, 2) + blitz::pow(k * dz - z0, 2));
  solution = where(
    r_global <= rr,
    1 + 0.25 * pow(1 + cos(pi * r_global / rr), 2),
    1.
  );

  run.advector(0) =  dt / dx;
  run.advector(1) =  0;
  run.advector(2) =  0;

  run.advance(nt);

  auto L2_error = sqrt(sum(pow(solution - run.advectee_global(), 2)));
  std::cout << "L2 error: " << L2_error << std::endl;
  if(L2_error > 20.3) throw std::runtime_error("L2 error greater than threshold (20.3)");
}

int main()
{
  {
    enum { opts = 0};
    enum { opts_iters = 1};
    test<opts, opts_iters>("upwind");
  }
}
