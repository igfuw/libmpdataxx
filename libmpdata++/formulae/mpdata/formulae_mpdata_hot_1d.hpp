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
      // 3rd order term (first term from eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998)
      template<opts_t opts, class arr_1d_t>
      inline auto HOT(
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::tot)>::type* = 0 
      ) -> decltype(0)
      { 
        return 0; //no Higher Order Terms for second accuracy 
      }

      template<opts_t opts, class arr_1d_t>
      inline auto HOT_helper(
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i
      ) return_macro(,
	(
          3 * GC(i+h) * abs(GC(i+h)) / ((formulae::G<opts>(G, i) + formulae::G<opts>(G, i+1)) / 2) 
          - 2 * pow(GC(i+h), 3) / pow((( formulae::G<opts>(G, i) + formulae::G<opts>(G, i+1)) / 2), 2)  
          - GC(i+h)
        ) / 3
      )

      template<opts_t opts, class arr_1d_t>
      inline auto HOT( // for positive sign signal (so no need for abs())
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::tot) && !opts::isset(opts, opts::abs) && !opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(,
        HOT_helper<opts>(GC, G, i) 
        * frac<opts>(
	    psi(i+2) - psi(i+1) - psi(i) + psi(i-1)
	    ,//-----------------------------------
	    psi(i+2) + psi(i+1) + psi(i) + psi(i-1)
	)
      )

      template<opts_t opts, class arr_1d_t>
      inline auto HOT( // for variable sign signal (hence the need for abs())
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::tot) && opts::isset(opts, opts::abs) && !opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
        HOT_helper<opts>(GC, G, i)
        * frac<opts>(
	    abs(psi(i+2)) - abs(psi(i+1)) - abs(psi(i)) + abs(psi(i-1))
	    ,//------------------------------------------------------- 
	    abs(psi(i+2)) + abs(psi(i+1)) + abs(psi(i)) + abs(psi(i-1))
	)
      )

      template<opts_t opts, class arr_1d_t>
      inline auto HOT( // for infinite gauge option (sum of psi -> sum of 1) 
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::tot) && opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs are mutually exclusive");
        ,
        HOT_helper<opts>(GC, G, i) 
        * (psi(i+2) - psi(i+1) - psi(i) + psi(i-1)) 
        / //--------------------------------------- 
        (1 + 1 + 1 + 1)
      )
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
