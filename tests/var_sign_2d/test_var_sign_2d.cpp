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

#include <libmpdata++/solvers/adv/donorcell_2d.hpp>
#include <libmpdata++/solvers/adv/mpdata_2d.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;
using real_t = float;

enum {x, y};
int nt = 10;

template <class solver_t, class slvs_t>
void add_solver(
  slvs_t &slvs, 
  const typename solver_t::real_t offset, 
  const std::string fname
) {
  typename solver_t::params_t p;

  // pre instantiation
  p.span = {24, 24};
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
  {
    std::ostringstream tmp;
    tmp << fname << "_offset=" << offset << "_%d_%d.svg";    
    p.gnuplot_output = tmp.str();
  }
  p.outvars = {{0, {.name = "psi", .unit = "1"}}};

  // instantiation
  slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(p));

  // post instantiation
  blitz::firstIndex i;
  blitz::secondIndex j;
  slvs.back().advectee() = offset + exp(
    -sqr(i+.5-p.span[x]/2.) / (2.*pow(p.span[x]/10, 2)) // TODO: assumes dx=dy=1
    -sqr(j+.5-p.span[y]/2.) / (2.*pow(p.span[y]/10, 2))
  );  
  slvs.back().advector(x) = .5; 
  slvs.back().advector(y) = .25;
}

int main() 
{
  boost::ptr_vector<concurr::any<real_t, 2>> slvs;

  for (auto &offset : std::vector<real_t>({0,-.5}))
  {
    add_solver<output::gnuplot<solvers::donorcell_2d<real_t>>>(
      slvs, offset, "donorcell"
    );
    add_solver<output::gnuplot<solvers::mpdata_2d<real_t>>>(
      slvs, offset, "mpdata_it=2"
    );
    add_solver<output::gnuplot<solvers::mpdata_2d<real_t, formulae::opts::pds>>>(
      slvs, offset, "mpdata-pds_it=2"
    );
  }

  for (auto &slv : slvs) slv.advance(nt);
}
