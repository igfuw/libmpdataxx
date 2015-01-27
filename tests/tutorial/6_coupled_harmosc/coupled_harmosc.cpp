/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a minimalistic model of a harmonic oscillator
 * (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
 */

#include "coupled_harmosc.hpp"
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx;

#include <boost/math/constants/constants.hpp>

// ------------------------------------------------------------
  template <class ct_params_t>
  struct coupled_harmosc_stats : coupled_harmosc<ct_params_t>
  {
    using parent_t = coupled_harmosc<ct_params_t>;
    using parent_t::parent_t;
   
    void hook_post_step()
    {
      parent_t::hook_post_step();
      this->mem->barrier();
      std::cerr << "AQQ:" << this->timestep << " rank=" << this->mem->rank() << std::endl;
    }
  };
// ------------------------------------------------------------


const int nt = 1400;

int main() 
{
//<listing-1>
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
    enum { rhs_scheme = 
      solvers::rhs_scheme_t::trapez };
    struct ix { enum {psi, phi}; };
  };
//</listing-1>
  using real_t = double;
//  using real_t = typename ct_params_t::real_t;
  
  using sim_t = output::gnuplot<
    coupled_harmosc_stats<ct_params_t>
  >;
  typename sim_t::rt_params_t p; 

//<listing-2>
  // run-time parameters
  using boost::math::constants::pi;
  p.dt = 1;
  p.omega = 2 * pi<real_t>() / p.dt / 400;
//</listing-2>
  p.grid_size = {1001};
  p.outfreq = 200; //10; 

  using ix = typename ct_params_t::ix;
  p.outvars = {
    {ix::psi, {.name = "psi", .unit = "1"}},
    {ix::phi, {.name = "phi", .unit = "1"}}
  };
  p.gnuplot_command = "plot";

  // instantiation
  concurr::threads<sim_t, bcond::cyclic, bcond::cyclic> run(p);

  // initial condition
  {
    blitz::firstIndex i;
    run.advectee(ix::psi) = where(        
      i<= 50 || i>= 150,                            // if
      0,                                            // then
      0.5 * (1 + cos(2 * pi<real_t>() * i / 100))   // else
    );
    run.advectee(ix::phi) = real_t(0);
  }
  //run.advector() = .5;
  run.advector() = 0.;

  //stats
 
/*
  int mult = 1;    //number of 1/4 oscillations
  float full = 100.;  //1/4 of full oscillation
*/
  // integration
  //run.advance(mult*full*p.dt);
  run.advance(100);

  //the sign doesn't matter cause we compare psi^2 and phi^2
//  decltype(run.advectee(ix::psi)) solution(run.advectee().extent());
//  if (mult % 2){
//  {
//    blitz::firstIndex i;
//    solution = where(        
//      i<= (50 /*+ 0.5*mult*full*/) || i>= (150 /*+ 0.5*mult*full*/),  // if
//      0,                                            // then
//      0.5 * (1 + cos(2 * pi<real_t>() * i / 100.))   // else
//    );
//   }}
//  else{
//  {
//    blitz::firstIndex i;
//    solution = where(        
//      i<= (50 /*+ 0.5*mult*full*/) || i>= (150 /*+ 0.5*mult*full*/),  // if
//      0,                                            // then
//      0.5 * (1 + cos(2 * pi<real_t>() * i / 100.))   // else
//    );
//   }}

//  decltype(run.advectee(ix::psi)) tmp(run.advectee().extent());
//  tmp = pow(pow(run.advectee(ix::psi), 2) +  pow(run.advectee(ix::phi), 2) - pow(solution, 2), 2);
//  std::cerr <<"rms = " << sqrt((sum(tmp) - tmp(tmp.extent(0)-1)) / (tmp.extent(0) - 1)) / mult / full / p.dt << std::endl;
  

  std::cerr<<std::fixed;
  std::cerr<<std::setprecision(20);
  std::cerr<<"min(psi) = "<< min(run.advectee(ix::psi))<<std::endl;
  std::cerr<<minIndex(run.advectee(ix::psi)) <<std::endl;
  std::cerr<<"max(psi) = "<< max(run.advectee(ix::psi))<<std::endl;
  std::cerr<<maxIndex(run.advectee(ix::psi)) <<std::endl;
  std::cerr<<"min(phi) = "<< min(run.advectee(ix::phi))<<std::endl;
  std::cerr<<minIndex(run.advectee(ix::phi)) <<std::endl;
  std::cerr<<"max(phi) = "<< max(run.advectee(ix::phi))<<std::endl;
  std::cerr<<maxIndex(run.advectee(ix::phi)) <<std::endl;

/*
  if ((mult % 4) == 0) {
    std::cerr<<maxIndex(run.advectee(ix::psi)) <<std::endl;
    std::cerr<<maxIndex(solution) <<std::endl;
  } 
  if ((mult % 4) == 1){
    std::cerr<<minIndex(run.advectee(ix::phi)) <<std::endl;
    std::cerr<<maxIndex(solution) <<std::endl;
  }
  if ((mult % 4) == 2) {
    std::cerr<<minIndex(run.advectee(ix::psi)) <<std::endl;
    std::cerr<<maxIndex(solution) <<std::endl;
  } 
  if ((mult % 4) == 3){
    std::cerr<<maxIndex(run.advectee(ix::phi)) <<std::endl;
    std::cerr<<maxIndex(solution) <<std::endl;
  }

  std::cerr <<"advectees = " << sum(pow(run.advectee(ix::psi), 2)) + sum(pow(run.advectee(ix::phi), 2))<<std::endl;
  std::cerr <<"solution = " << sum(pow(solution, 2))  << std::endl;
*/  
}
