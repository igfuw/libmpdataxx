/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>

namespace libmpdataxx
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      template<opts_t opts, class arr_1d_t>
      inline auto DFL(
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::dfl)>::type* = 0 
      ) -> decltype(0)
      { 
        return 0;  
      }

      template<opts_t opts, class arr_1d_t>
      inline auto DFL(
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::toa) && opts::isset(opts, opts::dfl)>::type* = 0 
      ) return_macro(,
        - GC(i+h) / (G(i+1) + G(i)) / 2 * (GC((i+1)+h) - GC(i-h))
      )
      // TODO: dfl && toa version of the above?
      
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
