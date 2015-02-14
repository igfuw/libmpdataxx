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

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
#include <set>

using namespace libmpdataxx;

using T = float;
int n = 501, nt = 1600;

T min[2] = {2, -1}, max[2] = {4, 1};

template <class T>
void setup(T &solver, int n) 
{
  blitz::firstIndex i;
  int width = 50, center = 100;
  solver.advectee(0) = where(i <= center-width/2 || i >= center+width/2, min[0], max[0]); 
  solver.advectee(1) = where(i <= center-width/2 || i >= center+width/2, min[1], max[1]); 
  solver.advector() = -.5; 
}

template <class T>
void setopts(T &p, const int nt, const std::string &fname, int n_iters)
{
  p.n_iters = n_iters;

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

template <opts::opts_t opt, class vec_t>
void add_solver(vec_t &slvs, const std::string &fname, int n_iters)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
    enum { opts = opt };
  };
  using output_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  typename output_t::rt_params_t p;
  setopts(p, nt, fname, n_iters);
  p.grid_size = {n};
  slvs.push_back(new concurr::threads<output_t, bcond::cyclic, bcond::cyclic>(p));
  setup(slvs.back(), n);
}

int main() 
{
  const int n_dims = 1;
  boost::ptr_vector<concurr::any<T, n_dims>> slvs, slvs_fct;

  add_solver<opts::abs>(slvs, "mpdata_iters=1", 1);
  add_solver<opts::abs>(slvs, "mpdata_iters=2", 2);
  add_solver<opts::abs>(slvs, "mpdata_iters=3", 3);

  add_solver<opts::abs | opts::fct>(slvs_fct, "mpdata_fct_iters=2", 2);
  add_solver<opts::abs | opts::fct>(slvs_fct, "mpdata_fct_iters=3", 3);

  // non-FCT solvers
  for (auto &slv : slvs) slv.advance(nt);

  // FCT solvers
  for (auto &slv : slvs_fct) 
  {
    slv.advance(nt);
    for (auto i : std::set<int>({0,1}))
    {
      T 
        mn = blitz::min(slv.advectee(i)),
	mx = blitz::min(slv.advectee(i));
      if (mn < min[i]) { std::cerr << mn << " < " << min[i] << std::endl; throw; }
      if (mx > max[i]) { std::cerr << mx << " > " << max[i] << std::endl; throw; }
    }
  }
}
