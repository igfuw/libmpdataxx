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
  int width = 50, center = 100;
  solver.state() = where(i <= center-width/2 || i >= center+width/2, 2, 4); 
  solver.courant() = .5; 
}

template <class T>
void setopts(T &p, const int nt, const std::string &fname)
{
  p.outfreq = nt; // siplays initial condition and the final state
  p.gnuplot_output = fname + ".svg";    
  p.outvars = {
    {0, {.name = "psi", .unit = "1"}}
  };
  p.gnuplot_command = "plot";
  p.gnuplot_with = "histeps";
  p.gnuplot_yrange = "[1.75:4.25]";
}

template <class solver_t, class vec_t>
void add_solver(vec_t &slvs, const std::string &fname)
{
  using output_t = output::gnuplot<solver_t>;
  typename output_t::params_t p;
  setopts(p, nt, fname);
  slvs.push_back(new concurr::threads<output_t, bcond::cyclic>(n, p));
  setup(slvs.back(), n);
}

int main() 
{
  boost::ptr_vector<concurr::any<real_t, 1>> slvs;

  add_solver<solvers::mpdata_1d<real_t, 1>>(slvs, "mpdata_iters=1");
  add_solver<solvers::mpdata_1d<real_t, 2>>(slvs, "mpdata_iters=2");
  add_solver<solvers::mpdata_1d<real_t, 3>>(slvs, "mpdata_iters=3");

  add_solver<solvers::mpdata_fct_1d<real_t, 2>>(slvs, "mpdata_fct_iters=2");
  add_solver<solvers::mpdata_fct_1d<real_t, 3>>(slvs, "mpdata_fct_iters=3");

  for (auto &slv : slvs) slv.advance(nt);
}
