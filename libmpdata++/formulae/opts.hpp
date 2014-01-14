/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace opts
    {
      using opts_t = unsigned long;

      // http://stackoverflow.com/questions/523724/c-c-check-if-one-bit-is-set-in-i-e-int-variable
      constexpr opts_t bit(const opts_t &x)
      {
	return opts_t(1) << x;
      }

      constexpr bool isset(const opts_t &x, const opts_t &y)
      {
        return 0 != (x & y);
      }

      enum
      {
        fct = opts::bit(0), // flux-corrected transport
        pds = opts::bit(1), // positive-definite signal (turns optimisations on)
        toa = opts::bit(2), // third-order accuracy terms
        eps = opts::bit(3), // use frac=nom/(den+eps) instead of frac=where(den!=0,nom/den,0) 
        npa = opts::bit(4), // use nprt=(x-abs(x))/2 instead of nprt=min(0,x), and analogous formulae for pprt
        iga = opts::bit(5), // infinite-gauge option
        nug = opts::bit(6), // non-unit G (default G = 1) - see Smolarkiewicz 2006 eq (25) and discussion below for info on G
        dfl = opts::bit(7)  // devergent flows
      };

    }; // namespace opts

    // nprt: implemented using min
    template<opts::opts_t opts, class arr_t>
    inline auto negpart(
      const arr_t &x,
      typename std::enable_if<!opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == false
    ) return_macro(,
      min(0, x)
    )

    // nprt: implemented using abs
    template<opts::opts_t opts, class arr_t>
    inline auto negpart(
      const arr_t &x,
      typename std::enable_if<opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == true
    ) return_macro(,
      (x - abs(x)) / 2
    )

    // pprt: implemented using max
    template<opts::opts_t opts, class arr_t>
    inline auto pospart(
      const arr_t &x,
      typename std::enable_if<!opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == false
    ) return_macro(,
      max(0, x)
    )

    // pprt: implemented using abx
    template<opts::opts_t opts, class arr_t>
    inline auto pospart(
      const arr_t &x,
      typename std::enable_if<opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == true
    ) return_macro(,
      (x + abs(x)) / 2
    )

    // 1D: G = const = 1
    template<opts::opts_t opts, class arr_t>
    inline typename arr_t::T_numtype G(
      const arr_t &G,
      const rng_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 2D: G = const = 1
    template<opts::opts_t opts, int d, class arr_t>
    inline typename arr_t::T_numtype G(
      const arr_t &G,
      const rng_t &, const rng_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 3D: G = const = 1
    template<opts::opts_t opts, int d, class arr_t>
    inline typename arr_t::T_numtype G(
      const arr_t &G,
      const rng_t &, const rng_t &, const rng_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 1D: G != const
    template<opts::opts_t opts, class arr_t> 
    inline auto G(
      const arr_t &G,
      const rng_t &i,
      typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
    ) return_macro(,
      G(i) + 0 // return_macro includes a call to blitz::safeToReturn() which expects an expression as an arg
    )

    // 2D: G != const
    template<opts::opts_t opts, int d, class arr_t> 
    inline auto G(
      const arr_t &G,
      const rng_t &i,
      const rng_t &j,
      typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
    ) return_macro(,
      G(idxperm::pi<d>(i, j)) + 0
    )

  }; // namespace formulae
}; // namespace libmpdataxx
