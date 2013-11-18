/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "gnuplot-iostream_1d/test_gnuplot-iostream_1d.cpp"
 * \image html "../../tests/gnuplot-iostream_1d/figure_iters=1.svg"
 * \image html "../../tests/gnuplot-iostream_1d/figure_iters=2.svg"
 * \image html "../../tests/gnuplot-iostream_1d/figure_iters=3.svg"
 */

#include <libmpdata++/solvers/adv/mpdata_1d.hpp>
#include <libmpdata++/solvers/adv/donorcell_1d.hpp>
#include <libmpdata++/solvers/adv/leapfrog_1d.hpp>
#include <libmpdata++/bcond/bcond.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

using real_t = float;
int n = 20, nt = 20;

template <class T>
void setup(T &solver, int n) 
{
  blitz::firstIndex i;
  solver.advectee() = exp(
    -sqr(.5+i-n/2.) / (2.*pow(n/10, 2)) // TODO: assumes dx=1
  );  
  solver.advector() = .5; 
}

template <class T>
void setopts(T &p, int nt, int n_iters)
{
  p.n_iters = n_iters;

  p.outfreq = nt / 10; 
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << n_iters << ".svg";
    p.gnuplot_output = tmp.str();    
  }
  p.outvars = {{0, {.name = "psi", .unit = "1"}}};
}

template <class slvs_t>
void add_solver(slvs_t &slvs, int it)
{
  using solver_t = output::gnuplot<solvers::mpdata_1d<real_t>>;
  typename solver_t::params_t p;
  setopts(p, nt, it);
  slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(n, p));
  setup(slvs.back(), n);
}

int main() 
{
  boost::ptr_vector<concurr::any<real_t, 1>> slvs;
  add_solver(slvs, 1);
  add_solver(slvs, 2);
  add_solver(slvs, 3);
  for (auto &slv : slvs) slv.advance(nt);
}
