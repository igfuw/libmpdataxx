/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "gnuplot-iostream_2d/test_gnuplot-iostream_2d.cpp"
 * \image html "../../tests/gnuplot-iostream_2d/figure.svg"
 */

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

#include <set>

enum {x, y};

int main() 
{
  using namespace libmpdataxx;

  int nt = 96;

  for (auto &n_iters : std::set<int>({1,2,4}))
  {
    // compile-time parameters
    struct ct_params_t : ct_params_default_t
    { 
      using real_t = float; 
      enum { n_dims = 2 };
      enum { n_eqs = 1 }; 
      enum { opts = 0 };
    };
    using solver_t = output::gnuplot<solvers::mpdata<ct_params_t>>;

    // run-time parameters
    solver_t::rt_params_t p;
    p.span = {24, 24};
    p.n_iters = n_iters;
    p.outfreq = nt;
    p.gnuplot_with = "lines";
    p.gnuplot_border = "4095";
    p.gnuplot_maxcolors = 42;
    p.gnuplot_zrange = "[-.666:1]";
    p.gnuplot_cbrange = "[-.025:1.025]";
    {
      std::ostringstream tmp;
      tmp << "figure_iters=" << n_iters << "_%s_%d.svg";
      p.gnuplot_output = tmp.str();    
    }
    p.outvars = {{0, {.name = "psi", .unit = "1"}}};

    // instantiation
    concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(p);

    // post-instantiation
    {
      blitz::firstIndex i;
      blitz::secondIndex j;
      slv.advectee() = exp(
	-sqr(.5+i-p.span[x]/2.) / (2.*pow(p.span[x]/10, 2)) // TODO: assumes dx=dy=1
	-sqr(.5+j-p.span[y]/2.) / (2.*pow(p.span[y]/10, 2)) 
      );  
      slv.advector(x) = .5; 
      slv.advector(y) = .25;
    }

    // calculations
    slv.advance(nt);
  } 
}
