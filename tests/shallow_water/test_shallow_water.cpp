/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <boost/math/constants/constants.hpp>

#include <advoocat/solvers/mpdata_2d.hpp>
#include <advoocat/bcond/cyclic_2d.hpp>
#include <advoocat/concurr/threads.hpp>
#include <advoocat/output/gnuplot.hpp>
using namespace advoocat; // TODO
#include "shallow_water.hpp"

enum {qx, qy, h};  // eqations
enum {x, y};       // dimensions

const int n_iters = 2;
const int nx = 32, ny = 32, nt = 100;

template <class T>
void setopts(T &p)
{
  p.dt = .05;
  p.dx = 1;
  p.dz = 1;

  p.outfreq = nt / 25;
  p.outvars = {
    {h,   {.name = "h",   .unit = "m"}}, 
  };
  p.gnuplot_with = "lines";
  p.gnuplot_output = "figure_%s_%d.svg";
  p.gnuplot_zrange = p.gnuplot_cbrange = "[.85:1.1]";
}

template <class T>
void setup(T &solver) 
{
  using boost::math::constants::pi;

  blitz::firstIndex i;
  blitz::secondIndex j;

  solver.state(h) = 1 - .1 * pow(
    sin(i * pi<typename T::real_t>() / nx) * 
    sin(j * pi<typename T::real_t>() / ny), 
    32
  );

  solver.state(qx) = 0;
  solver.state(qy) = 0;
}

int main() 
{
  using real_t = double;
  boost::ptr_vector<concurr::any<real_t, 2>> slvs;

  { 
    using solver_t = output::gnuplot<
      shallow_water<
	real_t, n_iters,
	qx, qy, h
      >
    >;
    solver_t::params_t p; 
    setopts(p);

    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
    setup(slvs.back());
  }

  for (auto &slv : slvs) slv.advance(nt); 
};
