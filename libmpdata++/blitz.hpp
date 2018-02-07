/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#if (!defined(__FAST_MATH__) || !defined(NDEBUG)) && !defined(BZ_DEBUG)
#  warning neither __FAST_MATH__ && NDEBUG nor BZ_DEBUG defined
#  warning   -ffast-math (Clang and GCC) and -DNDEBUG are recomended for release-mode builds
#  warning   -DBZ_DEBUG is recommended for debug-mode builds
#endif

#if defined(BZ_THREADSAFE)
#  error libmpdata++ uses blitz::neverDeleteData, please unset BZ_THREADSAFE
#endif

// force use of #pragma ivdep even if Blitz thinks the compiler does not support it
// (as of gcc 20140212, it gives an ICE: http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60198) - TODO: check in CMake
//#define BZ_USE_ALIGNMENT_PRAGMAS  

#include <blitz/tv2fastiter.h> // otherwise Clang fails in debug mode
#include <blitz/array.h>

  
#include <libmpdata++/kahan_reduction.hpp>

//////////////////////////////////////////////////////////
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/predef.h>

#if BOOST_COMP_GNUC || BOOST_COMP_CLANG
  #define forceinline_macro inline __attribute__((always_inline))
#else
  #define forceinline_macro inline
#endif

// C++11 auto return type macro, deprecated since C++14
#define return_macro(init,expr)          \
  -> decltype(blitz::safeToReturn(expr)) \
{                                        \
  init                                   \
  return safeToReturn(expr);             \
} 

namespace libmpdataxx
{
  template <int n_dims> using idx_t = blitz::RectDomain<n_dims>;
  using rng_t = blitz::Range;

  // non-int ix_t means either rng_t or idx_t
  template <class ix_t, class expr_t>
  forceinline_macro auto return_helper(const expr_t &expr, typename std::enable_if<!std::is_same<ix_t, int>::value>::type* = 0)
  {
    return blitz::safeToReturn(expr);
  }

  template <class ix_t, class expr_t>
  forceinline_macro auto return_helper(const expr_t &expr, typename std::enable_if<std::is_same<ix_t, int>::value>::type* = 0)
  {
    return expr;
  }

  // helper for getting the underlaying real_t from either arrays or scalar expressions
  // non-int ix_t means expr_t is a blitz array
  template <class ix_t, class expr_t>
  struct real_t_helper
  {
    using type = typename expr_t::T_numtype;
  };
  
  template <class expr_t>
  struct real_t_helper<int, expr_t>
  {
    using type = expr_t;
  };

  // Boost ptr_vector 
  template <class arr_t>
  struct arrvec_t : boost::ptr_vector<arr_t> 
  {
    using parent_t = boost::ptr_vector<arr_t>;

    const arr_t &operator[](const int i) const 
    {   
      return this->at(
	(i + this->size()) % this->size()
      );  
    }
    
    arr_t &operator[](const int i)
    {   
      return this->at(
        (i + this->size()) % this->size()
      );  
    }

    void push_back(arr_t *arr)
    {
      parent_t::push_back(arr);

#if !defined(NDEBUG)
      // filling the array with NaNs to ease debugging
      *arr = blitz::has_signalling_NaN(*arr->dataFirst())
	? blitz::signalling_NaN(*arr->dataFirst())
	: blitz::quiet_NaN(*arr->dataFirst());
#endif
    }
  };
} // namespace libmpdataxx

