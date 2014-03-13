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
using namespace libmpdataxx;

enum {x, y, z};
struct ct_params_t : ct_params_default_t
{
  using real_t = double;
//<listing-1>
  enum { n_dims = 3 };
  enum { n_eqs = 1 };
//</listing-1>
  enum { opts = formulae::opts::abs };
};

template<class T>
void setup(T &solver)
{
  const ct_params_t::real_t
    dt = 0.2,
    dx = 2.5,
    dy = 2.5,
    dz = 2.5,
    h = 4.,
    r = 7 * dx,
    x0 = (20 - 7 * pow(6, -0.5)) * dx,
    y0 = (20 - 7 * pow(6, -0.5)) * dy,
    z0 = (20 + 14 * pow(6, -0.5)) * dz;

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;

  // sphere shape
  decltype(solver.advectee()) tmp(solver.advectee().extent());
  tmp =   blitz::pow((i + 0.5) * dx - x0, 2)
        + blitz::pow((j + 0.5) * dx - y0, 2)
        + blitz::pow((k + 0.5) * dx - z0, 2);
  solver.advectee() = where(tmp - pow(r, 2) <= 0, h * (1 - blitz::sqrt(tmp) / r) , 0);

  const ct_params_t::real_t
    omega = 0.1,
    xc = 20 * dx,
    yc = 20 * dy,
    zc = 20 * dz;
  // constant angular velocity rotational field
  solver.advector(x) = (-omega * pow(2, -0.5) * ((j+.5) * dy - yc) + omega / 2 * ((k+.5) * dz - zc)) * dt / dx;
  solver.advector(y) = (omega * pow(2, -0.5) * ((i+.5) * dx - xc) - omega / 2 * ((k+.5) * dz - zc)) * dt / dy;
  solver.advector(z) = (-omega / 2 * ((i+.5) * dx - xc) + omega / 2 * ((j+.5) * dy - yc)) * dt / dz;
}

int main()
{
  int nt = 5 * 314;
//<listing-2>
  using solver_t = output::hdf5_xdmf<
    solvers::mpdata<ct_params_t>
  >;
//</listing-2>
  solver_t::rt_params_t p;

  // pre instantation
  p.n_iters = 4;
  p.span = {41, 41, 41};

  p.outfreq = nt;
  p.outvars[0].name = "psi";
  p.outdir = "test";

  // instantation
  concurr::threads<
  solver_t,
  bcond::cyclic, bcond::cyclic,
  bcond::cyclic, bcond::cyclic,
  bcond::cyclic, bcond::cyclic
  > slv(p);

  // post instantation
  setup(slv);

  // time stepping
  slv.advance(nt);
}
