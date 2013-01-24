/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "blitz.hpp"
#include "arakawa_c.hpp"
#include "phc.hpp"

namespace diagnose 
{
  // pressure as a function of "theta times dry air density"
  template <class arg_t, typename real_t>
  inline auto p(
    const arg_t &tht,
    const rng_t &i,
    const rng_t &j
  ) return_macro(
      phc::p_1000<real_t>() * real_t(pow(
      (tht(i,j) * si::kilograms / si::cubic_metres  * si::kelvins * phc::R_d<real_t>()) / phc::p_1000<real_t>(), 1 / (1 - phc::R_d_over_c_pd<real_t>())
    )) / si::pascals
  )
}; //diagnose

