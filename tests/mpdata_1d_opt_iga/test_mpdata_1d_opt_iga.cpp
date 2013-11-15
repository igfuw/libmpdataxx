/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief example showing how to use the iga option of mpdata 
 *
 * \include "mpdata_1d_opt_iga/test_mpdata_1d_opt_iga.cpp"
 * \image html "../../tests/mpdata_1d_opt_iga/figure_iters=1.svg" TODO
 * \image html "../../tests/mpdata_1d_opt_iga/figure_iters=2.svg" TODO
 * \image html "../../tests/mpdata_1d_opt_iga/figure_iters=3.svg" TODO
 */

#include <libmpdata++/solvers/adv/mpdata_1d.hpp>
#include <libmpdata++/solvers/adv/mpdata_fct_1d.hpp>
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
  solver.advectee(0) = where(i <= center-width/2 || i >= center+width/2, 1, 3); 
//  solver.advectee(1) = where(i <= center-width/2 || i >= center+width/2,  0, 2); 
  solver.advector() = .5; 
}

template <class T>
void setopts(T &p, const int nt, const std::string &fname)
{
  p.n_eqs = 1;
  p.outfreq = nt; // displays initial condition and the final state
  p.gnuplot_output = fname + ".svg";    
  p.outvars = {
    {0, {.name = "variable-sign signal", .unit = "1"}}//,
//    {1, {.name = "single-sign signal", .unit = "1"}}
  };
  p.gnuplot_command = "plot";
  p.gnuplot_with = "histeps";
  p.gnuplot_yrange = "[-2:3]";
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
  const int n_dims = 1;
  boost::ptr_vector<concurr::any<real_t, n_dims>> slvs;

  add_solver<solvers::mpdata_fct_1d<real_t, 2/*, formulae::opts::iga*/>>(slvs, "mpdata_iters=2");
//  add_solver<solvers::mpdata_1d<real_t, 2, formulae::opts::iga>>(slvs, "mpdata_iters=2_iga");

  for (auto &slv : slvs) slv.advance(nt);
}
