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
    // overloads of abs/min/max/where that pick out the correct version based on ix_t
    template<class ix_t, class arg_t>
    forceinline_macro auto abs(const arg_t &a, typename std::enable_if<std::is_same<ix_t, int>::value>::type* = 0)
    {
      return std::abs(a);
    }
    
    template<class ix_t, class arg_t>
    forceinline_macro auto abs(const arg_t &a, typename std::enable_if<std::is_same<ix_t, rng_t>::value>::type* = 0)
    {
      return blitz::abs(a);
    }

    template<class ix_t, class a_t, class b_t>
    forceinline_macro auto min(const a_t &a, const b_t &b, typename std::enable_if<std::is_same<ix_t, int>::value>::type* = 0)
    {
      using cm_t = typename std::common_type<a_t, b_t>::type;
      return std::min(static_cast<cm_t>(a), static_cast<cm_t>(b));
    }
    
    template<class ix_t, class a_t, class b_t>
    forceinline_macro auto min(const a_t &a, const b_t &b, typename std::enable_if<std::is_same<ix_t, rng_t>::value>::type* = 0)
    {
      return blitz::min(a, b);
    }

    template<class ix_t, class a_t, class b_t>
    forceinline_macro auto max(const a_t &a, const b_t &b, typename std::enable_if<std::is_same<ix_t, int>::value>::type* = 0)
    {
      using cm_t = typename std::common_type<a_t, b_t>::type;
      return std::max(static_cast<cm_t>(a), static_cast<cm_t>(b));
    }
    
    template<class ix_t, class a_t, class b_t>
    forceinline_macro auto max(const a_t &a, const b_t &b, typename std::enable_if<std::is_same<ix_t, rng_t>::value>::type* = 0)
    {
      return blitz::max(a, b);
    }
   
    // variadic max & min
    template<class ix_t, class a_t, class... b_ts>
    forceinline_macro auto max(const a_t &a, const b_ts & ... bs)
    {
      return max<ix_t>(a, max<ix_t>(bs...));
    }
    
    template<class ix_t, class a_t, class... b_ts>
    forceinline_macro auto min(const a_t &a, const b_ts & ... bs)
    {
      return min<ix_t>(a, min<ix_t>(bs...));
    }
      
    template<class ix_t, class arg_t>
    forceinline_macro auto where(bool c, const arg_t &a, const arg_t &b, typename std::enable_if<std::is_same<ix_t, int>::value>::type* = 0)   
    {
      return c ? a : b;
    }
    
    template<class ix_t, class c_t, class a_t, class b_t>
    forceinline_macro auto where(c_t c, const a_t &a, const b_t &b, typename std::enable_if<std::is_same<ix_t, rng_t>::value>::type* = 0)  
    {
      return blitz::where(c, a, b);
    }

    // nprt: implemented using min
    template<opts::opts_t opts, class ix_t, class arr_t>
    forceinline_macro auto negpart(
      const arr_t &x,
      typename std::enable_if<!opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == false
    )
    {
      return return_helper<ix_t>(min<ix_t>(0, x));
    }

    // nprt: implemented using abs
    template<opts::opts_t opts, class ix_t, class arr_t>
    forceinline_macro auto negpart(
      const arr_t &x,
      typename std::enable_if<opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == true
    )
    {
      return return_helper<ix_t>((x - abs<ix_t>(x)) / 2);
    }

    // pprt: implemented using max
    template<opts::opts_t opts, class ix_t, class arr_t>
    forceinline_macro auto pospart(
      const arr_t &x,
      typename std::enable_if<!opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == false
    )
    {
      return return_helper<ix_t>(max<ix_t>(0, x));
    }

    // pprt: implemented using abx
    template<opts::opts_t opts, class ix_t, class arr_t>
    forceinline_macro auto pospart(
      const arr_t &x,
      typename std::enable_if<opts::isset(opts, opts::npa)>::type* = 0 // enabled if npa == true
    )
    {
      return return_helper<ix_t>((x + abs<ix_t>(x)) / 2);
    }

    // 1D or ND: G = const = 1
    template<opts::opts_t opts, class arr_t, class ix_t>
    inline auto G(
      const arr_t &G,
      const ix_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 2D: G = const = 1
    template<opts::opts_t opts, int d, class arr_t, class ix_t>
    inline auto G(
      const arr_t &G,
      const ix_t &, const ix_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 3D: G = const = 1
    template<opts::opts_t opts, int d, class arr_t, class ix_t>
    inline auto G(
      const arr_t &G,
      const ix_t &, const ix_t &, const ix_t &,
      typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
    ) {
      return 1; 
    }

    // 1D on ND: G != const
    template<opts::opts_t opts, class arr_t, class ix_t>
    inline auto G(
      const arr_t &G,
      const ix_t &i,
      typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
    ) 
    {
      return return_helper<ix_t>(
        G(i) + 0 // return_macro includes a call to blitz::safeToReturn() which expects an expression as an arg
      );
    }
    
    // 2D: G != const
    template<opts::opts_t opts, int d, class arr_t, class ix_t> 
    inline auto G(
      const arr_t &G,
      const ix_t &i,
      const ix_t &j,
      typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
    )
    {
      return return_helper<ix_t>(
        G(idxperm::pi<d>(i, j)) + 0
      );
    }
    
    // 3D: G != const
    template<opts::opts_t opts, int d, class arr_t, class ix_t> 
    inline auto G(
      const arr_t &G,
      const ix_t &i,
      const ix_t &j,
      const ix_t &k,
      typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
    )
    {
      return return_helper<ix_t>(
        G(idxperm::pi<d>(i, j, k)) + 0
      );
    }
  } // namespace formulae
} // namespace libmpdataxx
