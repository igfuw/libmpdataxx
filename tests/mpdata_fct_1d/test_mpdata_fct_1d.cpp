/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief MPDATA-FCT example from @copybrief Smolarkiewicz_2006 (Figs 11-12)
 *
 * \include "mpdata_fct_1d/test_mpdata_fct_1d.cpp"
 * \image html "../../tests/mpdata_fct_1d/mpdata_iters=1.svg"
 * \image html "../../tests/mpdata_fct_1d/mpdata_iters=2.svg"
 * \image html "../../tests/mpdata_fct_1d/mpdata_fct_iters=2.svg"
 * \image html "../../tests/mpdata_fct_1d/mpdata_iters=3.svg"
 * \image html "../../tests/mpdata_fct_1d/mpdata_fct_iters=3.svg"
 */
// TODO: add sss runs to the list above

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
  solver.state(0) = where(i <= center-width/2 || i >= center+width/2, 2, 4); 
  solver.state(1) = where(i <= center-width/2 || i >= center+width/2, -1, 1); 
  solver.courant() = -.5; 
}

template <class T>
void setopts(T &p, const int nt, const std::string &fname)
{
  p.outfreq = nt; // diplays initial condition and the final state
  p.gnuplot_output = fname + ".svg";    
  p.outvars = {
    {0, {.name = "single-sign signal", .unit = "1"}},
    {1, {.name = "variable-sign signal", .unit = "1"}}
  };
  p.gnuplot_command = "plot";
  p.gnuplot_with = "histeps";
  p.gnuplot_yrange = "[-1.25:4.25]";
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

  const int n_eqs = 2;

//  add_solver<solvers::mpdata_1d<real_t, 1, n_eqs>>(slvs, "mpdata_iters=1");
//  add_solver<solvers::mpdata_1d<real_t, 2, n_eqs>>(slvs, "mpdata_iters=2");
//  add_solver<solvers::mpdata_1d<real_t, 3, n_eqs>>(slvs, "mpdata_iters=3");

  add_solver<solvers::mpdata_fct_1d<real_t, 2, n_eqs>>(slvs, "mpdata_fct_iters=2");
//  add_solver<solvers::mpdata_fct_1d<real_t, 3, n_eqs>>(slvs, "mpdata_fct_iters=3");

//  add_solver<solvers::mpdata_1d<real_t, 2, n_eqs, formulae::mpdata::sss>>(slvs, "mpdata_iters=2_sss");
//  add_solver<solvers::mpdata_fct_1d<real_t, 2, n_eqs, formulae::mpdata::sss>>(slvs, "mpdata_fct_iters=2_sss");

  for (auto &slv : slvs) slv.advance(nt);
}
