/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// advection (<> should be used instead of "" in normal usage) 
#include "advoocat/solvers/mpdata_2d.hpp"
#include "advoocat/bcond/cyclic_2d.hpp"
#include "advoocat/concurr/threads.hpp"
#include "advoocat/output/gnuplot.hpp"
using namespace advoocat; // TODO
#include "shallow_water.hpp"

enum {qx, qy, h};  // eqations
enum {x, y};       // dimensions

const int n_eqs = 3; // TODO not here!

using real_t = double;
const int n_iters = 2;

const int nx = 100, ny = 100, nt = 10;

template <class T>
void setopts(T &p)
{
  p.dt = .1;
  p.dHdx.resize(nx, ny);
  p.dHdx = 0;
  p.dHdx.resize(nx, ny);
  p.dHdy = 0;

  p.outvars = {
    {h,   {.name = "h",   .unit = "m"}}, 
  };
  p.gnuplot_output = "figure_%s_%d.svg";
}

int main() 
{

  boost::ptr_vector<concurr::any<real_t, 2>> slvs;

  { 
    using solver_t = output::gnuplot<
      shallow_water<
	solvers::mpdata_2d<real_t, n_iters, n_eqs>,
	qx, qy, h
      >
    >;
    solver_t::params_t p; 
    setopts(p);

    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }

  for (auto &slv : slvs)
  {
    // TODO: initial condition

    // integration
    slv.advance(nt); 
  }
};
