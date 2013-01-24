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
  namespace detail 
  {
    // pressure as a function of "theta times dry air density"
    template <typename real_t>
    inline quantity<si::pressure, real_t> p(
      const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th 
    )
    {  
      return phc::p_1000<real_t>() * real_t(pow(  // real_t() needed due to a bug in Boost.units [1]
        rhod_th * phc::R_d<real_t>() / phc::p_1000<real_t>(), 
        1 / (1 - phc::R_d_over_c_pd<real_t>())
      ));
    }

    template <typename real_t>
    inline real_t p_dimles(const real_t rhod_th)
    {
      return p(rhod_th * si::kelvins * si::kilograms / si::cubic_metres) / si::pascals;
    }
  }

  float p(float rhod_th) { return detail::p_dimles(rhod_th); }
  double p(double rhod_th) { return detail::p_dimles(rhod_th); }
  long double p(long double rhod_th) { return detail::p_dimles(rhod_th); }

  BZ_DECLARE_FUNCTION(p) // works with Blitz > 0.10 [2]

  // [1] https://svn.boost.org/trac/boost/ticket/6957
  // [2] http://sourceforge.net/mailarchive/message.php?msg_id=30394213

}; //diagnose

