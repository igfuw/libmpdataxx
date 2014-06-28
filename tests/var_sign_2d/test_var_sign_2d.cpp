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

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;
using T = float;

enum {x, y};
int nt = 10;

template <class solver_t, class slvs_t>
void add_solver(
  slvs_t &slvs, 
  const typename solver_t::real_t offset, 
  const std::string fname
) {
  typename solver_t::rt_params_t p;

  // pre instantiation
  p.grid_size = {24, 24};
  p.outfreq = nt; 
  p.gnuplot_with = "lines";
  p.gnuplot_border = "4095";
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
  slvs.push_back(new concurr::threads<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  >(p));

  // post instantiation
  blitz::firstIndex i;
  blitz::secondIndex j;
  slvs.back().advectee() = offset + exp(
    -sqr(i-(p.grid_size[x]-1)/2.) / (2.*pow((p.grid_size[x]-1)/10, 2)) // TODO: assumes dx=dy=1
    -sqr(j-(p.grid_size[y]-1)/2.) / (2.*pow((p.grid_size[y]-1)/10, 2))
  );  
  slvs.back().advector(x) = .5; 
  slvs.back().advector(y) = .25;
}

int main() 
{
  boost::ptr_vector<concurr::any<T, 2>> slvs;

  for (auto &offset : std::vector<T>({0,-.5}))
  {
    {
      struct ct_params_t : ct_params_default_t
      {
        using real_t = T;
        enum { n_dims = 2 };
        enum { n_eqns = 1 };
        enum { opts = opts::abs };
      };
      add_solver<output::gnuplot<solvers::mpdata<ct_params_t>>>(
        slvs, offset, "mpdata-abs_it=2"
      );
    }
    {
      struct ct_params_t : ct_params_default_t
      {
        using real_t = T;
        enum { n_dims = 2 };
        enum { n_eqns = 1 };
        enum { opts = 0 };
      };
      add_solver<output::gnuplot<solvers::mpdata<ct_params_t>>>(
        slvs, offset, "mpdata_it=2"
      );
    }
  }

  for (auto &slv : slvs) slv.advance(nt);
}
