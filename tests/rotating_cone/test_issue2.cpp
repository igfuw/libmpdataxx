/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "rotating_cone/test_rotating_cone.cpp"
 * \image html "../../tests/rotating_cone/figure.svg"
 */

#include <cmath>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
using namespace libmpdataxx;

template <typename ct_params_t>
struct solver_t : solvers::mpdata<ct_params_t>
{
  using parent_t = solvers::mpdata<ct_params_t>;

  // inheriting ctors
  using parent_t::parent_t;    

  // adding a hook
  void hook_post_step()
  {
    parent_t::hook_post_step();

    this->mem->barrier();
    if (this->mem->rank() != 0) return;
    if (min(this->mem->advectee()) < 0) 
    {
      auto ix = minIndex(this->mem->advectee());
      std::cerr << "  min(" << ix << ") = " << this->mem->advectee()(ix) << std::endl;
      throw std::runtime_error("values below zero!");
    }
  }
};

template <int opts_arg>
void test(const std::string filename, const int &nx, const int &ny)
{
  enum {x, y};
  struct ct_params_t : ct_params_default_t
  {
    using real_t = float;
    enum { n_dims = 2 };
    enum { n_eqs = 1 };
    enum { opts = opts_arg };
  };

  typename ct_params_t::real_t 
    dt = .1,
    dx = 1,
    dy = 1,
    omega = .1,
    h = 4., // TODO: other name!
    h0 = 0;

  int nt = 628 * 6;

  using sim_t = solver_t<ct_params_t>;
  typename sim_t::rt_params_t p;

  // pre instantiation
  p.n_iters = 2;
  p.span = {nx, ny};

  // instantiation
  concurr::threads<sim_t, bcond::open, bcond::open, bcond::open, bcond::open> run(p); 

  // post instantiation
  {
    typename ct_params_t::real_t
      r = 15. * dx,
      x0 = 75 * dx,
      y0 = 50 * dy,
      xc = .5 * p.span[x] * dx,
      yc = .5 * p.span[y] * dy;

    blitz::firstIndex i;
    blitz::secondIndex j;

    // cone shape
    decltype(run.advectee()) tmp(run.advectee().extent());
    tmp = blitz::pow((i+.5) * dx - x0, 2) + blitz::pow((j+.5) * dy - y0, 2);
    run.advectee() = h0 + where(tmp - pow(r, 2) <= 0, h * blitz::sqr(1 - tmp / pow(r, 2)), 0.);

    // constant angular velocity rotational field
    run.advector(x) = -omega * ((j+.5) * dy - yc) * dt / dx;
    run.advector(y) =  omega * ((i+.5) * dx - xc) * dt / dy;
    // TODO: an assert confirming that the above did what it should have done
    //       (in context of the advector()'s use of blitz::Array::reindex())
  }

  // time stepping
  run.advance(nt);
}

int main(int argc, char **argv)
{
  int 
    nx = atoi(argv[1]),
    ny = atoi(argv[2]);

  std::cerr << "nx=" << nx << " ny=" << ny << std::endl;

/*
  {
    enum { opts = 0 };
    std::cerr << "starting basic run" << std::endl;
    test<opts>("basic", nx, ny);
  }
*/
  {
    enum { opts = formulae::opts::fct };
    std::cerr << "starting fct run" << std::endl;
    test<opts>("fct", nx, ny);
  }
/*
  {
    enum { opts = formulae::opts::iga | formulae::opts::fct};
    std::cerr << "starting iga_fct run" << std::endl;
    test<opts>("iga_fct", nx, ny);
  }
  {
    enum { opts = formulae::opts::iga | formulae::opts::tot | formulae::opts::fct};
    std::cerr << "starting iga_fct_tot run" << std::endl;
    test<opts>("iga_fct_tot", nx, ny);
  }
*/
}
