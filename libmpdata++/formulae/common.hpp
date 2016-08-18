/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/idxperm.hpp>
#include <libmpdata++/opts.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
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
    inline int G(
      const arr_t &G,
      const rng_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 2D: G = const = 1
    template<opts::opts_t opts, int d, class arr_t>
    inline int G(
      const arr_t &G,
      const rng_t &, const rng_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 3D: G = const = 1
    template<opts::opts_t opts, int d, class arr_t>
    inline int G(
      const arr_t &G,
      const rng_t &, const rng_t &, const rng_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }
    
    // ND: G = const = 1
    template<opts::opts_t opts, class arr_t, class ijk_t>
    inline int G(
      const arr_t &G,
      const ijk_t &,
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
    
    // 3D: G != const
    template<opts::opts_t opts, int d, class arr_t> 
    inline auto G(
      const arr_t &G,
      const rng_t &i,
      const rng_t &j,
      const rng_t &k,
      typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
    ) return_macro(,
      G(idxperm::pi<d>(i, j, k)) + 0
    )
    
    // ND: G != const
    template<opts::opts_t opts, class arr_t, class ijk_t> 
    inline auto G(
      const arr_t &G,
      const ijk_t &ijk,
      typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
    ) return_macro(,
      G(ijk) + 0
    )
  } // namespace formulae
} // namespace libmpdataxx
