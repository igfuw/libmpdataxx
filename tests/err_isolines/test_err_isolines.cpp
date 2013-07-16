/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "err_isolines/test_err_isolines.cpp"
 * \image html "../../tests/err_isolines/figure.svg" 
 */

#include <fstream>
#include <list>

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#include <boost/math/constants/constants.hpp>

#include <libmpdata++/solvers/adv/donorcell_1d.hpp>
#include <libmpdata++/solvers/adv/mpdata_1d.hpp>
#include <libmpdata++/solvers/adv/mpdata_fct_1d.hpp>
#include <libmpdata++/bcond/bcond.hpp>
#include <libmpdata++/concurr/threads.hpp>


// making things simpler (yet less elegant)
using namespace libmpdataxx;
using boost::math::constants::pi;
blitz::firstIndex i;
boost::ptr_map<std::string, std::ofstream> outfiles;
using real_t = long double; // this is a good test to show differences betwee float and double!!!


// helper function template to ease adding the solvers to the pointer map
template <class solver_t, class vec_t>
void add_solver(vec_t &slvs, const std::string &key, const int nx)
{
  using params_t = typename solver_t::params_t;

  boost::assign::ptr_map_insert<
    concurr::threads<solver_t, bcond::cyclic> // map element type
  >(slvs)(
    key, // map key
    nx, params_t() // concurr's ctor args
  );

  boost::assign::ptr_map_insert(outfiles)(key);
  if (!outfiles[key].is_open())
  {
    char bits[3];
    sprintf(bits, "%03lu", 8 * sizeof(real_t));
    outfiles[key].open("err_mpdata_" + key + "_" + bits + "_bits.txt", std::ios::trunc); 
  }
}


// all the test logic
int main() 
{
  // gauss shape functor definition 
  struct gauss_t
  {
    // member fields
    real_t A0, A, sgma, x0;

    // call operator
    real_t operator()(real_t x) const 
    { 
      return A0 + A * exp( real_t(-.5) * pow(x - x0, 2) / pow(sgma, 2));
    }
    
    // Blitz magick
    BZ_DECLARE_FUNCTOR(gauss_t)
  };

  // simulation parameters
  const real_t 
    t_max    = 1., // "arbitrarily"
    dx_max   = 1.,
    x_max    = 44. * dx_max, // see not about compact support in asserts below
    sgma     = 1.5 * dx_max, 
    velocity = dx_max / t_max, // "solution advects over the one grid increment for r=8"
    x0       = .5 * x_max, 
    A0       = -.5,
    A        = 1. / sgma / sqrt(2 * pi<real_t>());

  const std::list<real_t> 
    courants({ .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95}),
    dxs({dx_max, dx_max/2, dx_max/4, dx_max/8, dx_max/16, dx_max/32, dx_max/64, dx_max/128 });

  // gauss shape functor instantiation
  gauss_t gauss({.A0 = A0, .A = A, .sgma = sgma, .x0 = x0});

  // looping over different grid increments
  for (auto &dx : dxs) 
  {
    std::cerr << "dx = " << dx << std::endl;
    
    // looping over different Courant numbers
    for (auto &cour : courants)
    { 
      std::cerr << "  C = " << cour << std::endl;

      const real_t dt = cour * dx / velocity;
      const int 
        n_dims = 1,
        n_eqs  = 1,
        nx = round(x_max/dx),
        nt = round(t_max/dt);

      boost::ptr_map<std::string, concurr::any<real_t, n_dims>> slvs;

      // silly loop order, but it helped to catch a major bug!

      // donor-cell
      add_solver<solvers::donorcell_1d<real_t, n_eqs>>(slvs, "iters=1", nx);

      // MPDATA
      add_solver<solvers::mpdata_1d<real_t, 2, n_eqs>>(slvs, "iters=2", nx);
//      add_solver<solvers::mpdata_1d<real_t, 2, n_eqs, formulae::mpdata::toa>>(slvs, "iters=2_toa", nx);
      add_solver<solvers::mpdata_1d<real_t, 3, n_eqs>>(slvs, "iters=3", nx);
      add_solver<solvers::mpdata_1d<real_t, 3, n_eqs, formulae::mpdata::toa | formulae::mpdata::iga>>(slvs, "iters=3_toa", nx);

      // MPDATA-FCT
//      add_solver<solvers::mpdata_fct_1d<real_t, 2, n_eqs>>(slvs, "iters=2_fct", nx);
//      add_solver<solvers::mpdata_fct_1d<real_t, 2, n_eqs, formulae::mpdata::toa>>(slvs, "iters=2_fct_toa", nx);
//      add_solver<solvers::mpdata_fct_1d<real_t, 3, n_eqs>>(slvs, "iters=3_fct", nx);
//      add_solver<solvers::mpdata_fct_1d<real_t, 3, n_eqs, formulae::mpdata::toa>>(slvs, "iters=3_fct_toa", nx);

      // calculating the analytical solution
      typename solvers::donorcell_1d<real_t, n_eqs>::arr_t exact(nx);
      exact = gauss((i+.5)*dx - velocity * dt * nt);

      // looping over solvers
      for (auto keyval : slvs) 
      {
        auto &key = keyval.first;
        auto &slv = *keyval.second;

        std::cerr << "    solver = " << key << std::endl; 

        // setting the solver up
	slv.courant() = cour; 
        slv.state() = gauss((i+.5)*dx);
   
        // running the solver
	slv.advance(nt);

        // asserting that periodic boundries did not affect the result
        // and that the chosen domain length is enough to have compact support up to machine precision
        assert(exact(0) == slv.state()(0));
        assert(exact(nx-1) == slv.state()(nx-1));

        // calculating the deviation from analytical solution
        real_t err = sqrt(sum(pow(slv.state() - exact, 2)) / nx) / (nt * dt);

        outfiles[key] << dx << "\t" << cour << "\t" << err << std::endl;
      }
    }
  }
}
