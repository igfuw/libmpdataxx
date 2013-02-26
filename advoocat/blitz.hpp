/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/
/** @mainpage
 *  TODO: table of examples with columns:
 *  - adv. algorithm
 *  - Jacobian
 *  - number of dimensions
 *  - equation set
 *
 *  blitz-based
 * 
 *  suggested compiler options (by compiler): -march=native, -Ofast, -DNDEBUG, -lblitz (opt), -DBZDEBUG (opt), -std=c++11
 */

// code licensed under the terms of GNU GPL v3
// copyright holder: University of Warsaw

#pragma once

#if defined(_OPENMP) || defined(_REENTRANT)
#  define BZ_THREADSAFE
#endif
#include <blitz/array.h>

// local definition of a Kahan's sum reduction
// (http://en.wikipedia.org/wiki/Kahan_summation_algorithm)
namespace blitz
{
  template<typename P_sourcetype, typename P_resulttype = BZ_SUMTYPE(P_sourcetype)>
  class ReduceKahanSum 
  {
    public:

    typedef P_sourcetype T_sourcetype;
    typedef P_resulttype T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceKahanSum() { } 

    bool operator()(const T_sourcetype& x,const int=0) const 
    { 
      T_resulttype y, t;
// TODO: check which level
#pragma GCC push_options
#pragma GCC optimize ("O0")
      y = x - c_;
      t = sum_ + y;
      c_ = (t - sum_) - y;
      sum_ = t;
#pragma GCC pop_options
      return true;
    }   

    T_resulttype result(const int) const { return sum_; }

    void reset() const 
    { 
      sum_ = c_ = zero(T_resulttype()); 
    }
 
    static const char* name() { return "sum"; }
 
    protected:

    mutable T_resulttype sum_, c_;
  };
  BZ_DECL_ARRAY_PARTIAL_REDUCE(kahan_sum, ReduceKahanSum)
  BZ_DECL_ARRAY_FULL_REDUCE(kahan_sum, ReduceKahanSum)
}
  
#include <boost/ptr_container/ptr_vector.hpp>

// C++11 auto return type macro
#define return_macro(init,expr)          \
  -> decltype(blitz::safeToReturn(expr)) \
{                                        \
  init                                   \
  return safeToReturn(expr);             \
} 

namespace advoocat
{
  template <int n_dims> using idx_t = blitz::RectDomain<n_dims>;
  using rng_t = blitz::Range;

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

}; // namespace advoocat
