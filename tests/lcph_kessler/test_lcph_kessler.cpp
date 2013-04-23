#include <advoocat/bcond/cyclic_2d.hpp>
#include <advoocat/concurr/threads.hpp>
#include <advoocat/output/gnuplot.hpp>

#include "cloud.hpp"

// for initial condition
#include <libcloudph++/common/phc_hydrostatic.hpp>
#include <libcloudph++/common/phc_theta.hpp>

using namespace advoocat;

enum {rhod_th_ix, rhod_rv_ix, rhod_rc_ix, rhod_rr_ix }; // variables
enum {x, z}; // dimensions
using real_t = double;

// 8th ICMW case 1 by Wojciech Grabowski)
#include "icmw8_case1.hpp"


// simulation and output parameters
template <class T>
void setopts(T &p, int nt, int n_iters)
{
  //p.outfreq = nt/10; // TODO
  //p.gnuplot_zrange = p.gnuplot_cbrange = "[.5:2.5]";
  p.gnuplot_view = "map";
  {
    std::ostringstream tmp;
    tmp << "output/figure_iters=" << n_iters << "_%s_%d.svg";
    p.gnuplot_output = tmp.str();
  }
  p.outvars = 
  {
    {rhod_rv_ix, {.name = "\\rho_v", .unit = "kg/m^{-3}"}},
    {rhod_rc_ix, {.name = "\\rho_c", .unit = "kg/m^{-3}"}}
  };
}


int main()
{
  int nx = 32, nz = 32, nt = 20;
  const int n_iters = 2;

  // helper type to shorten the code below
  using solver_t = output::gnuplot<cloud<real_t, n_iters, 
    rhod_th_ix, 
    rhod_rv_ix,
    rhod_rc_ix,
    rhod_rr_ix
  >>;

  // instantiation of structure containing simulation parameters
  solver_t::params_t p;

  // output and simulation parameters
  setopts(p, nt, n_iters);

  // solver instantiation
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(nx, nz, p);

  // initial condition
  icmw8_case1::setup(slv);

  // timestepping
  slv.advance(nt);
}
