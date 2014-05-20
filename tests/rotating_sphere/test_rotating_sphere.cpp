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
//</listing-1>
  enum { n_eqns = 1 };
  enum { opts = opts::abs };
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
  tmp =   blitz::pow(i * dx - x0, 2)
        + blitz::pow(j * dx - y0, 2)
        + blitz::pow(k * dx - z0, 2);
  solver.advectee() = where(tmp - pow(r, 2) <= 0, h * (1 - blitz::sqrt(tmp) / r) , 0);

  const ct_params_t::real_t
    omega = 0.1,
    xc = 20 * dx,
    yc = 20 * dy,
    zc = 20 * dz;
  // constant angular velocity rotational field
  solver.advector(x) = omega / sqrt(3) * (-(j * dy - yc) + (k * dz - zc)) * dt / dx;
  solver.advector(y) = omega / sqrt(3) * ( (i * dx - xc) - (k * dz - zc)) * dt / dy;
  solver.advector(z) = omega / sqrt(3) * (-(i * dx - xc) + (j * dy - yc)) * dt / dz;
}

int main()
{
  int nt = 5 * 314;
  using slv_t = solvers::mpdata<ct_params_t>;
//<listing-2>
  using slv_out_t = output::hdf5_xdmf<slv_t>;
//</listing-2>
  slv_out_t::rt_params_t p;

  // pre instantation
  p.n_iters = 4;
  p.grid_size = {41, 41, 41};

  p.outfreq = nt;
  p.outvars[0].name = "psi";
//<listing-3>
  p.outdir = "rotating_sphere_3d";
//</listing-3>

  // instantation
  concurr::threads<
  slv_out_t,
  bcond::open, bcond::open,
  bcond::open, bcond::open,
  bcond::open, bcond::open
  > slv(p);

  // post instantation
  setup(slv);

  // time stepping
  slv.advance(nt);
}
