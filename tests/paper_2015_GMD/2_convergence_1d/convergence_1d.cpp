/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * convergence test 
 */

#include <fstream>
#include <list>

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#include <boost/math/constants/constants.hpp>

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/cxx11_thread.hpp>

// making things simpler (yet less elegant)
using namespace libmpdataxx;
using boost::math::constants::pi;
blitz::firstIndex i;
boost::ptr_map<std::string, std::ofstream> outfiles;
using T = double; // with long double this is a good test to show differences between float and double!!!

// helper function template to ease adding the solvers to the pointer map
template <opts::opts_t opt, class vec_t>
void add_solver(vec_t &slvs, const std::string &key, const int nx, const int n_iters)
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

  p.n_iters = n_iters;
  p.grid_size = {nx};

  boost::assign::ptr_map_insert<
    concurr::cxx11_thread<solver_t, bcond::open, bcond::open> // map element type
  >(slvs)(
    key,  // map key
    p     // concurr's ctor args
  );

  boost::assign::ptr_map_insert(outfiles)(key);
  if (!outfiles[key].is_open())
  {
    outfiles[key].open("err_mpdata_" + key + ".txt", std::ios::trunc); 
  }
}

// integrated gauss shape
struct gauss_int_t 
{
  //member fields
  T A0, A, sgma, x0, dx;

  //integral from a to b of gauss function
  T operator()(T a) const
  {
    return A0 * dx + A * sgma * sqrt(pi<T>()/2) * (erf( (a+dx/2 - x0) / sqrt(2) / sgma ) - erf((a-dx/2 - x0) / sqrt(2) / sgma));
  } 

  // Blitz magick
  BZ_DECLARE_FUNCTOR(gauss_int_t)
};

// all the test logic
int main() 
{
  // simulation parameters
  const T 
    t_max    = 1.,                // "arbitrarily"
    dx_max   = 1.,
    x_max    = 10 * 44. * dx_max, // see note about compact support in asserts below
    sgma     = 1.5 * dx_max, 
    velocity = dx_max / t_max,    // "solution advects over the one grid increment for r=8"
    x0       = .5 * x_max, 
    A0       = 0,
    A        = 1. / sgma / sqrt(2 * pi<T>());

  const std::list<T> 
    courants({ .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95}),
    dxs({dx_max, dx_max/2, dx_max/4, dx_max/8, dx_max/16, dx_max/32, dx_max/64, dx_max/128 });

  // looping over different grid increments
  for (auto &dx : dxs) 
  {
    std::cerr << "dx = " << dx << std::endl;

    // gauss shape functor instantiation
    gauss_int_t gauss_int({.A0 = A0, .A = A, .sgma = sgma, .x0 = x0, .dx = dx});

    // looping over different Courant numbers
    for (auto &cour : courants)
    { 
      std::cerr << "  C = " << cour << std::endl;

      const T dt = cour * dx / velocity;
      const int 
        n_dims = 1,
        nx = round(x_max/dx),
        nt = round(t_max/dt);

      boost::ptr_map<std::string, concurr::any<T, n_dims>> slvs;

      // silly loop order, but it helped to catch a major bug!

      // donor-cell
      add_solver<0>(slvs, "iters=1", nx, 1);

      // MPDATA
      add_solver<0>(slvs, "iters=2", nx, 2);
      //add_solver<opts::tot>(slvs, "iters=2_tot", nx, 2);
      add_solver<0>(slvs, "iters=3", nx, 3);
      add_solver<opts::tot>(slvs, "iters=3_tot", nx, 3);
      add_solver<opts::iga>(slvs, "iters=i", nx, 2);
      //add_solver<opts::iga | opts::tot>(slvs, "iters=i_tot", nx, 2);

      // MPDATA-FCT
      add_solver<opts::fct>(slvs, "iters=2_fct", nx, 2);
      //add_solver<opts::fct | opts::tot>(slvs, "iters=2_fct_tot", nx, 2);
      //add_solver<opts::fct>(slvs, "iters=3_fct", nx, 3);
      //add_solver<opts::fct | opts::tot>(slvs, "iters=3_fct_tot", nx, 3);
      add_solver<opts::fct | opts::iga>(slvs, "iters=i_fct", nx, 2);
      add_solver<opts::fct | opts::iga | opts::tot>(slvs, "iters=i_fct_tot", nx, 2);

      // calculating the analytical solution
      decltype(slvs.end()->second->advectee()) exact(nx);
      exact = gauss_int(i*dx - velocity * dt * nt) / dx;

      // looping over solvers
      for (auto keyval : slvs) 
      {
        auto &key = keyval.first;
        auto &slv = *keyval.second;

        std::cerr << "    solver = " << key << std::endl; 

        // setting the solver up
	slv.advector() = cour; 
        slv.advectee() = gauss_int(i*dx) / dx;
   
        // running the solver
	slv.advance(nt);

        // asserting that boundary conditions do not affect the result
        // and that the chosen domain length is enough to have compact support up to machine precision
        // exact(0) === 0; slv.advectee(0) === 0;
        assert(exact(0) == slv.advectee()(0));
        assert(exact(nx-1) == slv.advectee()(nx-1));

        // calculating the deviation from analytical solution
        T err = sqrt(sum(pow(slv.advectee() - exact, 2)) / nx) / (nt * dt);

        outfiles[key] << std::scientific << std::setprecision(4) <<std::endl;
        outfiles[key] << dx << "\t" << cour << "\t" << err << std::endl;
      }
    }
  }
}
