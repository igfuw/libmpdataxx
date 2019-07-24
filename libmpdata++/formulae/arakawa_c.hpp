/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/blitz.hpp>

namespace libmpdataxx
{
  namespace arakawa_c
  {
    namespace
    {
      struct hlf_t {} h; 

      rng_t operator+(
        const rng_t &i, const hlf_t &
      ) { 
        return i; 
      } 

      rng_t operator-(
        const rng_t &i, const hlf_t &
      ) { 
        return i-1; 
      }
      
      int operator+(
        const int i, const hlf_t &
      ) { 
        return i; 
      } 

      int operator-(
        const int i, const hlf_t &
      ) { 
        return i-1; 
      }

      template<class n_t>
      rng_t operator^(
        const rng_t &r, const n_t &n
      ) { 
        return rng_t(
          (r - n).first(), 
          (r + n).last()
        ); 
      } 
    };
  } // namespace arakawa_c
} // namespace libmpdataxx
