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

#include <blitz/array.h>

template <int n_dims> using idx_t = blitz::RectDomain<n_dims>;
using rng_t = blitz::Range;

// C++11 auto return type macro
#define return_macro(expr) \
  -> decltype(blitz::safeToReturn(expr)) \
{ return safeToReturn(expr); } 

// Boost ptr_vector 
#include <boost/ptr_container/ptr_vector.hpp>
template <class arr_t>
struct arrvec_t : boost::ptr_vector<arr_t> 
{
  using parent = boost::ptr_vector<arr_t>;

  const arr_t &operator[](const int i) const 
  {   
    return this->at(
      (i + this->size()) % this->size()
    );  
  }

  void push_back(arr_t *arr)
  {
    parent::push_back(arr);

#if !defined(NDEBUG)
    // filling the array with NaNs to ease debugging
    *arr = blitz::has_signalling_NaN(*arr->dataFirst())
      ? blitz::signalling_NaN(*arr->dataFirst())
      : blitz::quiet_NaN(*arr->dataFirst());
#endif
  }
};

