/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "mpdata_1d_opt_npa/test_mpdata_1d_opt_eps.cpp"
 * \image html "../../tests/mpdata_1d_opt_npa/figure_iters=3.svg" TODO
 */

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

// TODO: make a common file with the setopts and setup from below?

using namespace libmpdataxx;

using T = float;
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
void setopts(T &p, const int nt, const std::string &fname, int n_iters)
{
  p.n_iters = n_iters;

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

template <formulae::opts::opts_t opt, class vec_t>
void add_solver(vec_t &slvs, const std::string &fname, int n_iters)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { n_dims = 1 };
    enum { n_eqs = 2 };
    enum { opts = opt };
  };
  using output_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  typename output_t::rt_params_t p;
  setopts(p, nt, fname, n_iters);
  p.span = {n};
  slvs.push_back(new concurr::threads<output_t, bcond::cyclic, bcond::cyclic>(p));
  setup(slvs.back(), n);
}

int main() 
{
  const int n_dims = 1;
  boost::ptr_vector<concurr::any<T, n_dims>> slvs;

  add_solver<formulae::opts::abs>(slvs, "mpdata_iters=2", 2);
  add_solver<formulae::opts::abs | formulae::opts::npa>(slvs, "mpdata_iters=2_npa", 2);
  add_solver<formulae::opts::abs>(slvs, "mpdata_iters=3", 3);
  add_solver<formulae::opts::abs | formulae::opts::npa>(slvs, "mpdata_iters=3_npa", 3);

  for (auto &slv : slvs) slv.advance(nt);
}
