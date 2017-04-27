/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the G factor

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // interpolation of G to (i+1/2) - general case
      template<opts_t opts, class arr_1d_t>
      inline auto G_bar_x( 
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0
      ) return_macro(,
        (
          formulae::G<opts>(G, i+1) + formulae::G<opts>(G, i)
        ) / 2
      )

      // interpolation of G to (i+1/2) - constant G version
      template<opts_t opts, class arr_1d_t>
      inline typename arr_1d_t::T_numtype G_bar_x(
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0
      ) {
        return 1;
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
