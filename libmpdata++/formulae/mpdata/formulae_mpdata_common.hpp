/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/idxperm.hpp>
#include <libmpdata++/arakawa_c.hpp>

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
        sss = detail::bit(0), // single-sign signal
        toa = detail::bit(1)//, // third-order accuracy // TODO (not code to handle it yet)
// TODO        eps = detail::bit(2)  // use frac=nom/(den+eps) instead of frac=where(den!=0,nom/den,0) // TODO! (and value of eps)
// TODO: gauge option (adding a large constant)
// TODO: negpart form option
      };

      constexpr bool opt_set(const opts_t &x, const opts_t &y) 
      { 
	return 0 != (x & y); 
      }

      const int /*halo = 1,*/ n_tlev = 2;

      constexpr const int halo(const opts_t &opts) 
      {
        return opt_set(opts, toa) ? 2 : 1; 
      }

      template<class nom_t, class den_t>
      inline auto frac(
        const nom_t &nom, const den_t &den
      ) return_macro(,
        where(den > 0, nom / den, 0)
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
