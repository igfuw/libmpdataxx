/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief MPDATA-FCT example from @copybrief Smolarkiewicz_2006 (Figs 11-12)
 *
 * \include "gnuplot-iostream_1d/test_gnuplot-iostream_1d.cpp"
 * \image html "../../tests/gnuplot-iostream_1d/figure_iters=1.svg"
 * \image html "../../tests/gnuplot-iostream_1d/figure_iters=2.svg"
 * \image html "../../tests/gnuplot-iostream_1d/figure_iters=3.svg"
 */

#include <libmpdata++/solvers/mpdata_fct_1d.hpp>
#include <libmpdata++/bcond/bcond.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

using real_t = float;
int n = 500, nt = 1600;

template <class T>
void setup(T &solver, int n) 
{
  blitz::firstIndex i;
  solver.state() = where(i < 100-12 || i > 100+12, 2, 4); 
  solver.courant() = .5; 
}

template <class T>
void setopts(T &p, int nt, int n_iters)
{
  p.outfreq = nt; // siplays initial condition and the final state
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << n_iters << ".svg";
    p.gnuplot_output = tmp.str();    
  }
  p.outvars = {
    {0, {.name = "psi", .unit = "1"}}
  };
  p.gnuplot_command = "plot";
}

template <int it, class slvs_t>
void add_solver(slvs_t &slvs)
{
  using solver_t = output::gnuplot<
    solvers::mpdata_fct_1d<real_t, it>
  >;
  typename solver_t::params_t p;
  setopts(p, nt, it);
  slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(n, p));
  setup(slvs.back(), n);
}

int main() 
{
  boost::ptr_vector<concurr::any<real_t, 1>> slvs;
  add_solver<1>(slvs);
  add_solver<2>(slvs);
  add_solver<3>(slvs);
  for (auto &slv : slvs) slv.advance(nt);
}
