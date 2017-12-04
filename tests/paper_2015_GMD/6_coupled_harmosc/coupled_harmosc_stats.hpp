/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * outputing min/max of the advected signal, calculating rms error
 */

#pragma once
#include <boost/math/constants/constants.hpp>
#include <iomanip>
#include <iostream>
#include <fstream>

template <class parent_t>
struct stats : public parent_t
{
  using parent_t::parent_t;
  using real_t = typename parent_t::real_t;
  using ix = typename parent_t::ix;

  std::ofstream ofs;                
 
  void hook_post_step()
  {
    parent_t::hook_post_step();
    this->mem->barrier();
    if (this->rank != 0) return;

    if (this->timestep % 100 == 0)
    {
      if (!ofs.is_open())
        ofs.open("stats.txt", std::ofstream::out);
      // outputing signal extrema ...
      ofs<< std::fixed << std::setprecision(8)<<std::endl;
      ofs<<"timestep = "<<this->timestep << std::endl;
      ofs<<"min(psi) = "<<std::setw(12) << min(this->mem->advectee(ix::psi))<<" @ index "<< minIndex(this->mem->advectee(ix::psi)) << std::endl;
      ofs<<"max(psi) = "<<std::setw(12) << max(this->mem->advectee(ix::psi))<<" @ index "<< maxIndex(this->mem->advectee(ix::psi)) << std::endl;
      ofs<<"min(phi) = "<<std::setw(12) << min(this->mem->advectee(ix::phi))<<" @ index "<< minIndex(this->mem->advectee(ix::phi)) << std::endl;
      ofs<<"max(phi) = "<<std::setw(12) << max(this->mem->advectee(ix::phi))<<" @ index "<< maxIndex(this->mem->advectee(ix::phi)) << std::endl;

      decltype(this->mem->advectee(ix::psi)) solution(this->mem->advectee().lbound(), this->mem->advectee().extent());
  
      assert(max(this->mem->advector() == min(this->mem->advector())));
      real_t vel = this->mem->advector()(this->mem->grid_size[0].first());

      real_t t = this->timestep * this->dt;
      using boost::math::constants::pi;

      // solution @ 1/4T, 1/2T, 3/4T and T
      // the sign doesn't matter cause we are going to compare to psi^2 and phi^2 
      // only one variable, cause the other should be equal zero at those times
      if (this->timestep % 2)
      {
        blitz::firstIndex i;
        solution = where(
          i<= (50 + vel * t) || i>= (150 + vel * t),    // if
          0,                                            // then
          0.5 * (1 + cos(2 * pi<real_t>() * i / 100.))  // else
        );                   
      }
      else
      {
        blitz::firstIndex i;
        solution = where(        
          i<= (50 + vel * t) || i>= (150 + vel * t),    // if
          0,                                            // then
          0.5 * (1 + cos(2 * pi<real_t>() * i / 100.))  // else
        );                  
      }

      // ... and rms error (normalised by the domain size and simulation time)
      decltype(this->mem->advectee(ix::psi)) tmp(this->mem->advectee().extent());
      tmp = pow(pow(this->mem->advectee(ix::psi), 2) +  pow(this->mem->advectee(ix::phi), 2) - pow(solution, 2), 2);
      ofs << "rms = " << std::scientific << sqrt((sum(tmp) - tmp(tmp.extent(0)-1)) / (tmp.extent(0) - 1)) / this->timestep << std::endl;
    }
  }
};
 
