/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief example showing how to use the pds option of mpdata intended for
 *   optimising calculations on fields of positive sign
 *
 * \include "mpdata_1d_opt_pds/test_mpdata_1d_opt_pds.cpp"
 * \image html "../../tests/mpdata_1d_opt_pds/figure_iters=1.svg" TODO
 * \image html "../../tests/mpdata_1d_opt_pds/figure_iters=2.svg" TODO
 * \image html "../../tests/mpdata_1d_opt_pds/figure_iters=3.svg" TODO
 */

#include <libmpdata++/solvers/adv/mpdata_1d.hpp>
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
  solver.advectee(0) = where(i <= center-width/2 || i >= center+width/2, -1, 1); 
  solver.advectee(1) = where(i <= center-width/2 || i >= center+width/2,  2, 4); 
  solver.advector() = .5; 
}

template <class T>
void setopts(T &p, const int nt, const std::string &fname)
{
  p.outfreq = nt; // displays initial condition and the final state
  p.gnuplot_output = fname + ".svg";    
  p.outvars = {
    {0, {.name = "variable-sign signal", .unit = "1"}},
    {1, {.name = "single-sign signal", .unit = "1"}}
  };
  p.gnuplot_command = "plot";
  p.gnuplot_with = "histeps";
  p.gnuplot_yrange = "[-2:5]";
}

template <class solver_t, class vec_t>
void add_solver(vec_t &slvs, const std::string &fname)
{
  using output_t = output::gnuplot<solver_t>;
  typename output_t::params_t p;
  setopts(p, nt, fname);
  p.span = {n};
  slvs.push_back(new concurr::threads<output_t, bcond::cyclic>(p));
  setup(slvs.back(), n);
}

int main() 
{
  const int n_dims = 1;
  boost::ptr_vector<concurr::any<real_t, n_dims>> slvs;

  add_solver<solvers::mpdata_1d<real_t, /* n_eqs = */ 2>>(slvs, "mpdata_iters=2");
  add_solver<solvers::mpdata_1d<real_t, /* n_eqs = */ 2, formulae::opts::pds>>(slvs, "mpdata_iters=2_pds");

  // TODO: test if pds gives any speed-up with single-sign field

  for (auto &slv : slvs) slv.advance(nt);
}
