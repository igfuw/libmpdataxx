/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "blitz.hpp"
#include "pi.hpp"
#include "arakawa_c.hpp"

namespace courant 
{
  template<int d, class arr_2d_t, typename real_t>
  void intrp(
    const arr_2d_t &C,
    const arr_2d_t &psi,
    const rng_t &im,
    const rng_t &j, 
    const real_t &dt, 
    const real_t &dx 
  ) {
    C(pi<d>(im+h,j)) = dt / dx * .5 * (psi(pi<d>(im,j)) + psi(pi<d>(im+1,j)));
  } 
}; //courant
