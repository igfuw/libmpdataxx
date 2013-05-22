/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "gnuplot-iostream_2d/test_gnuplot-iostream_2d.cpp"
 * \image html "../../tests/gnuplot-iostream_2d/figure.svg"
 */

#include <libmpdata++/solvers/mpdata_2d.hpp>
#include <libmpdata++/solvers/donorcell_2d.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

enum {x, y};

template <class T>
void setup(T &solver, int n[2]) 
{
  blitz::firstIndex i;
  blitz::secondIndex j;
  solver.state() = exp(
    -sqr(.5+i-n[x]/2.) / (2.*pow(n[x]/10, 2)) // TODO: assumes dx=dy=1
    -sqr(.5+j-n[y]/2.) / (2.*pow(n[y]/10, 2)) 
  );  
  solver.courant(x) = .5; 
  solver.courant(y) = .25;
}

template <class T>
void setopts(T &p, int nt, int n_iters)
{
  p.outfreq = nt;
  p.gnuplot_with = "lines";
  p.gnuplot_border = "4095";
  p.gnuplot_maxcolors = 42;
  p.gnuplot_zrange = "[-.666:1]";
  p.gnuplot_cbrange = "[-.025:1.025]";
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << n_iters << "_%s_%d.svg";
    p.gnuplot_output = tmp.str();    
  }
  p.outvars = {{0, {.name = "psi", .unit = "1"}}};
}

int main() 
{
  using namespace libmpdataxx;

  int n[] = {24, 24}, nt = 96;

  {
    using solver_t = output::gnuplot<solvers::donorcell_2d<float>>;
    solver_t::params_t p;
    setopts(p, nt, 1);
    concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p);

    setup(slv, n);
    slv.advance(nt);
  } 
  {
    const int it = 2;
    using solver_t = output::gnuplot<solvers::mpdata_2d<float, it>>;
    solver_t::params_t p;
    setopts(p, nt, it);
    concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p); 

    setup(slv, n); 
    slv.advance(nt);
  } 
  {
    const int it = 4;
    using solver_t = output::gnuplot<solvers::mpdata_2d<float, it>>;
    solver_t::params_t p;
    setopts(p, nt, it);
    concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p); 

    setup(slv, n); 
    slv.advance(nt); 
  }
}
