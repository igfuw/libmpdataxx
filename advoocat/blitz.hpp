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
 */

// code licensed under the terms of GNU GPL v3
// copyright holder: University of Warsaw

#pragma once

#include <blitz/array.h>

template <typename real_t>
using arr_1d_t = blitz::Array<real_t, 1>;

template <typename real_t>
using arr_2d_t = blitz::Array<real_t, 2>;

using rng_t = blitz::Range;
using idx_1d_t = blitz::RectDomain<1>;
using idx_2d_t = blitz::RectDomain<2>;

#define return_macro(expr) \
  -> decltype(blitz::safeToReturn(expr)) \
{ return safeToReturn(expr); } 

#include <boost/ptr_container/ptr_vector.hpp>
template <class arr_t>
struct arrvec_t : boost::ptr_vector<arr_t> {
  const arr_t &operator[](const int i) const {   
    return this->at(
      (i + this->size()) % this->size()
    );  
  }
};

