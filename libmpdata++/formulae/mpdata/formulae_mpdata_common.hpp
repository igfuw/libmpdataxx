/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/idxperm.hpp>
#include <libmpdata++/arakawa_c.hpp>

#include <boost/preprocessor/control/if.hpp>

namespace libmpdataxx
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      using namespace arakawa_c;
      using idxperm::pi;

      using opts_t = unsigned long;

      namespace detail
      {
        // http://stackoverflow.com/questions/523724/c-c-check-if-one-bit-is-set-in-i-e-int-variable
        constexpr opts_t bit(const opts_t &x) 
        { 
          return opts_t(1) << x;  
        }
      };

      enum 
      {
        pds = detail::bit(0), // positive-definite signal (turns optimisations on)
        toa = detail::bit(1), // third-order accuracy terms
        eps = detail::bit(2), // use frac=nom/(den+eps) instead of frac=where(den!=0,nom/den,0) 
        npa = detail::bit(3), // use nprt=(x-abs(x))/2 instead of nprt=min(0,x), and analogous formulae for pprt
        iga = detail::bit(4)  // infinite-gauge option
      };

      constexpr bool opt_set(const opts_t &x, const opts_t &y) 
      { 
	return 0 != (x & y); 
      }

      const int n_tlev = 2;

      constexpr const int halo(const opts_t &opts) 
      {
        return opt_set(opts, toa) ? 2 : 1; 
      }

      // frac: implemented using blitz::where()
      template<opts_t opts, class nom_t, class den_t>
      inline auto frac(
        const nom_t &nom, 
        const den_t &den,
        typename std::enable_if<!opt_set(opts, eps)>::type* = 0 // enabled if eps == false
      ) return_macro(,
        where(den > 0, nom / den, 0)
      )

      // frac: implemented as suggested in MPDATA papers
      template<opts_t opts, class nom_t, class den_t>
      inline auto frac(
        const nom_t &nom, 
        const den_t &den,
        typename std::enable_if<opt_set(opts, eps)>::type* = 0 // enabled if eps == true
      ) return_macro(,
        nom / (den + blitz::tiny(typename den_t::T_numtype(0))) // note: for negative signal eps -> -eps
      )

      // nprt: implemented using min
      template<opts_t opts, class arr_t> 
      inline auto negpart(
        const arr_t &x,
        typename std::enable_if<!opt_set(opts, npa)>::type* = 0 // enabled if npa == false
      ) return_macro(,
        min(0, x)
      )

      // nprt: implemented using abs
      template<opts_t opts, class arr_t> 
      inline auto negpart(
        const arr_t &x,
        typename std::enable_if<opt_set(opts, npa)>::type* = 0 // enabled if npa == true
      ) return_macro(,
        (x - abs(x)) / 2
      )

      // pprt: implemented using max
      template<opts_t opts, class arr_t> 
      inline auto pospart(
        const arr_t &x,
        typename std::enable_if<!opt_set(opts, npa)>::type* = 0 // enabled if npa == false
      ) return_macro(,
        max(0, x)
      )

      // pprt: implemented using abx
      template<opts_t opts, class arr_t> 
      inline auto pospart(
        const arr_t &x,
        typename std::enable_if<opt_set(opts, npa)>::type* = 0 // enabled if npa == true
      ) return_macro(,
        (x + abs(x)) / 2
      )
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
