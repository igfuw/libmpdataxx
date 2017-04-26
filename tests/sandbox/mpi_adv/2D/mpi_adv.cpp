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
    enum { n_dims = 2 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };

  int nx = 257;
  int ny = 65;

  T pi = boost::math::constants::pi<T>();
  T dx = 4.0 / (nx - 1);
  T dy = 1.0 / (ny - 1);
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
  p.grid_size = {nx, ny};

  p.outfreq = nt / 8; 
  p.outvars[0].name = "psi";
  p.outdir = "out";

  // instantiation
  concurr::threads<
    slv_out_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p); 

  T rr = 0.5;
  T x0 = 0.5;
  T y0 = 0.5;

  decltype(run.advectee()) r(run.advectee().extent());
  r.reindexSelf(run.advectee().lbound());

  blitz::firstIndex i;
  blitz::secondIndex j;

  r = sqrt(blitz::pow(i * dx - x0, 2) + blitz::pow(j * dy - y0, 2));
  run.advectee() = where(
    r <= rr,
    1 + 0.25 * pow(1 + cos(pi * r / rr), 2),
    1.
  );

  decltype(run.advectee()) solution(run.advectee().shape());
  solution.reindexSelf(run.advectee().lbound());

  T disp = std::fmod(time, 4.0);
  
  r = sqrt(blitz::pow(i * dx - (x0 + disp), 2) + blitz::pow(j * dy - y0, 2));
  solution = where(
    r <= rr,
    1 + 0.25 * pow(1 + cos(pi * r / rr), 2),
    1.
  );

  run.advector(0) =  dt / dx;
  run.advector(1) =  0;

  run.advance(nt);

  std::cout << "L2 error: " << sqrt(sum(pow(solution - run.advectee(), 2))) << std::endl;
}

int main()
{
  {
    enum { opts = 0};
    enum { opts_iters = 1};
    test<opts, opts_iters>("upwind");
  }
}
