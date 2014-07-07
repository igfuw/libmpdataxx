/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "mpdata_1d_opt_nug/test_mpdata_1d_opt_nug.cpp"
 * \image html "../../tests/mpdata_1d_opt_nug/figure_iters=3.svg" TODO
 */

#include <libmpdata++/solvers/mpdata.hpp> 
#include <libmpdata++/concurr/threads.hpp>
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream.h> 

using namespace libmpdataxx;

using T = float;
using arr_t = blitz::Array<T, 1>;

int n = 15;
arr_t G(n), phi(n); 
arr_t G_c(n+1);

template <class T>
void setopts(T &p, const std::string &sfx)
{
  p.n_iters =   2; //number of iterations
}

template <opts::opts_t opt, class vec_t>
void add_solver(vec_t &slvs, const std::string &sfx)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { n_dims = 1 };
    enum { n_eqns = 1 };
    enum { opts = opt };
  };
  using solver_t = solvers::mpdata<ct_params_t>;
  typename solver_t::rt_params_t p;
  setopts(p, sfx);
  p.grid_size = {n};
  slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(p));
}

template <class T0, class T1, class T2>
void plot(T0 &gp, T1 &slvs, T2 &G)
{
  gp << "plot" 
    << " '-' with steps title 'mxr = psi / G'"
    << ",'-' with steps title 'mxr = psi'\n";
  {
    arr_t tmp(slvs.at(0).advectee() / G);
    gp.send(tmp);
  }
  gp.send(slvs.at(1).advectee());
  gp << "plot" 
    << " '-' with steps title 'rho = psi'"
  << ",'-' with steps title 'rho = psi * G'\n";
  gp.send(slvs.at(0).advectee());
  gp.send(arr_t(slvs.at(1).advectee() * G));
}

int main() 
{
  const int n_dims = 1;
  boost::ptr_vector<concurr::any<T, n_dims>> slvs;

  // Starting point: a non-divergent G \times C
  T GC = .25;

  phi =   1,  5, 10, 10,  10,  5,  1,  1,   1,   1,   1,  1,  1,  1,  1;
  G   =   1,  1,  1,  1,   1,  1,  1,  1,  .5,  .5,  .5, .5, .5,  1,  1;
  G_c = 1,  1,  1,  1,   1,  1,  1,  1, .75,  .5,  .5, .5, .5, .75,  1,  1;      

  add_solver<0>(slvs, "");

  // advecting density with Courant number
  slvs.back().advectee() = G * phi;
  slvs.back().advector() = GC / G_c; 

  add_solver<opts::nug>(slvs, "_nug");
  // advecting mixing ratio with C times rho
  slvs.back().g_factor() = G;
  slvs.back().advectee() = phi;
  slvs.back().advector() = GC; 

  if (sum(G * phi) != sum(slvs.at(0).advectee())) 
    throw std::runtime_error("t=0:  G*phi != psi");
  if (sum(G * phi) != sum(slvs.at(1).advectee() * slvs.at(1).g_factor())) 
    throw std::runtime_error("t=0:  G*phi != psi*G");

  Gnuplot gp;
  gp << "set term svg size 1200, 500\n";
  gp << "set output 'output.svg'\n";
  gp << "set yrange [0:12]\n";
  gp << "set multiplot layout 2,3 columnsfirst\n";
  gp << "set grid\n";

  plot(gp, slvs, G);
  for (auto &slv : slvs) slv.advance(20);
  plot(gp, slvs, G);
  for (auto &slv : slvs) slv.advance(30);
  plot(gp, slvs, G);

  T eps = 1e-6;
  // have to use std::abs(), as abs() alone returns an integer!
  if (std::abs(sum(G * phi) - sum(slvs.at(0).advectee())) > eps) 
    throw std::runtime_error("t>0:  G*phi - psi > eps");
  if (std::abs(sum(G * phi) - sum(slvs.at(1).advectee() * slvs.at(1).g_factor())) > eps) 
    throw std::runtime_error("t>0:  G*phi - psi*G > eps");
}
