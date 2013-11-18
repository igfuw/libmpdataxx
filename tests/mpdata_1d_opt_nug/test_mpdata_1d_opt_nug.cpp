/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "mpdata_1d_opt_nug/test_mpdata_1d_opt_nug.cpp"
 * \image html "../../tests/mpdata_1d_opt_nug/figure_iters=3.svg" TODO
 */

#include <libmpdata++/solvers/adv/donorcell_1d.hpp> // TODO: this test is aimed at testing MPDATA
#include <libmpdata++/bcond/bcond.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

using real_t = float;
using arr_t = blitz::Array<real_t, 1>;

int n = 12, nt = 12; //n / C;
arr_t G(n), phi(n); 
arr_t G_c(n+1);

template <class T>
void setopts(T &p, const int nt, const std::string &fname)
{
  p.outfreq = nt; // displays initial condition and the final state
  p.gnuplot_output = fname + ".svg";    
  p.outvars = {{0, {.name = "psi", .unit = "1"}}};
  p.gnuplot_command = "plot";
  p.gnuplot_with = "steps";
  p.gnuplot_lt = "3";
  p.gnuplot_grid = false;
  //p.gnuplot_yrange = "[-2:5]";
}

template <class solver_t, class vec_t>
void add_solver(vec_t &slvs, const std::string &fname)
{
  using output_t = output::gnuplot<solver_t>;
  typename output_t::params_t p;
  setopts(p, nt, fname);
  slvs.push_back(new concurr::threads<output_t, bcond::cyclic>(n, p));
}

int main() 
{
  const int n_dims = 1;
  boost::ptr_vector<concurr::any<real_t, n_dims>> slvs;

  // Starting point: a non-divergeng G \times C
  real_t GC = .25;

  G   =   1,  1,  1,  1,  .5,  .5,  .5, .5,  1,  1,  1,  1;
  G_c = 1,  1,  1,  1,  .75, .5, .5,  .5, .75, 1,  1,  1,  1;      
  phi =   1, 10, 10,  1,   1,   1,   1,  1,  1,  1,  1,  1;

  add_solver<solvers::donorcell_1d<real_t>>(slvs, "mpdata_iters=1");
  // advecting density with Courant number
  slvs.back().g_factor() = 1;
  slvs.back().advectee() = G * phi;
  slvs.back().advector() = GC / G_c; 

  add_solver<solvers::donorcell_1d<real_t, formulae::opts::nug>>(slvs, "mpdata_iters=1_nug");
  // advecting mixing ratio with C times rho
  slvs.back().g_factor() = G;
  slvs.back().advectee() = phi;
  slvs.back().advector() = GC; 

  if (sum(G * phi) != sum(slvs.at(0).advectee())) throw;
  if (sum(G * phi) != sum(slvs.at(1).advectee() * slvs.at(1).g_factor())) throw;

  for (auto &slv : slvs) slv.advance(nt);

  real_t eps = 1e-22;
  if (abs(sum(G * phi) - sum(slvs.at(0).advectee())) > eps) throw;
  if (abs(sum(G * phi) - sum(slvs.at(1).advectee() * slvs.at(1).g_factor())) > eps) throw;
}
