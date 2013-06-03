/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief test for advection of variable sign field 
 *
 * \include "var_sign_2d/test_var_sign_2d.cpp"
 * \image html "../../tests/var_sign_2d/figure.svg"
 */

#include <libmpdata++/solvers/donorcell_2d.hpp>
#include <libmpdata++/solvers/mpdata_2d.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

enum {x, y};

template <class T>
void setup(T &solver, int n[2], typename T::real_t offset) 
{
  blitz::firstIndex i;
  blitz::secondIndex j;
  solver.state() = offset + exp(
    -sqr(i+.5-n[x]/2.) / (2.*pow(n[x]/10, 2)) // TODO: assumes dx=dy=1
    -sqr(j+.5-n[y]/2.) / (2.*pow(n[y]/10, 2))
  );  
  solver.courant(x) = .5; 
  solver.courant(y) = .25;
}

template <class T>
void setopts(T &p, const int nt, const float offset, const std::string &fname)
{
  p.outfreq = nt; 
  p.gnuplot_with = "lines";
  p.gnuplot_border = "4095";
  p.gnuplot_maxcolors = 42; 
  p.gnuplot_zrange = "[-.666:1]";
  {
    std::ostringstream tmp;
    tmp << "[" << (offset -.025) << ":" << (offset +1.025) << "]";
    p.gnuplot_cbrange = tmp.str();
  }
  p.gnuplot_output = fname + "_%d_%d.svg";    
  p.outvars = {{0, {.name = "psi", .unit = "1"}}};
}

int main() 
{
  using namespace libmpdataxx;

  int n[] = {24, 24}, nt = 10;

  for (auto &offset : std::vector<float>({0,-.5}))
  {
    {
      using solver_t = output::gnuplot<solvers::donorcell_2d<float>>;
      solver_t::params_t p;
      setopts(p, nt, offset, "donorcell");
      concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p);
      setup(slv, n, offset);
      slv.advance(nt);
    } 
    {
      using solver_t = output::gnuplot<solvers::mpdata_2d<float, 2>>;
      solver_t::params_t p;
      setopts(p, nt, offset, "mpdata_it=2");
      concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p); 
      setup(slv, n, offset); 
      slv.advance(nt);
    } 
    {
      const int n_eqs = 1;
      using solver_t = output::gnuplot<solvers::mpdata_2d<float, 2, n_eqs, formulae::mpdata::sss>>;
      solver_t::params_t p;
      setopts(p, nt, offset, "mpdata-sss_it=2");
      concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p); 
      setup(slv, n, offset); 
      slv.advance(nt); 
    }
  }
}
