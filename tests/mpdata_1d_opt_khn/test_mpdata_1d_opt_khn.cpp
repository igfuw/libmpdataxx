/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "mpdata_1d_opt_eps/test_mpdata_1d_opt_khn.cpp"
 * \image html "../../tests/mpdata_1d_opt_eps/figure_iters=3.svg" TODO
 */

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

// TODO: make a common file with the setopts and setup from below?

using namespace libmpdataxx;

using T = float;
int n = 21, nt = 20;

template <class slv_t>
void setup(slv_t &solver, int n) 
{
  blitz::firstIndex i;
  int width = 2, center = 3;
  solver.advectee() = where(i <= center-width/2 || i >= center+width/2, -40, 40) * blitz::tiny(T(0)); 
  solver.advector() = .5; 
}

template <class T>
void setopts(T &p, const int nt, const std::string &fname, int n_iters, const std::string &title)
{
  p.n_iters = n_iters;

  p.outfreq = nt; // displays initial condition and the final state
  p.gnuplot_output = fname + ".svg";    
  p.gnuplot_command = "plot";
  p.gnuplot_size = "ratio 3";
  p.gnuplot_title = title;
}

template <opts::opts_t opt, class vec_t>
void add_solver(vec_t &slvs, const std::string &fname, int n_iters, const std::string &title)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { n_dims = 1 };
    enum { n_eqns = 1 };
    enum { opts = opt };
  };
  using output_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  typename output_t::rt_params_t p;
  setopts(p, nt, fname, n_iters, title);
  p.grid_size = {n};
  slvs.push_back(new concurr::threads<output_t, bcond::cyclic, bcond::cyclic>(p));
  setup(slvs.back(), n);
}

int main() 
{
  const int n_dims = 1;
  boost::ptr_vector<concurr::any<T, n_dims>> slvs;

  add_solver<opts::abs>(slvs, "kahan_off", 2, "opts::abs");
  add_solver<opts::abs | opts::khn>(slvs, "kahan_on", 2, "opts::abs | opts::khn");

  for (auto &slv : slvs) slv.advance(nt);
}
