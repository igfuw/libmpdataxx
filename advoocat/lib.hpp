/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/
/** @mainpage
 * List of examples (subset of list of tests):
 * - test_gnuplot-iostream.cpp
 * - test_var_sign_2d.cpp
 */

// code licensed under the terms of GNU GPL v3
// copyright holder: University of Warsaw

typedef double real_t;

#include <blitz/array.h>
using arr_1d_t = blitz::Array<real_t, 1>;
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

struct hlf_t {} h;

inline rng_t operator+(
  const rng_t &i, const hlf_t &
) { 
  return i; 
} 

inline rng_t operator-(
  const rng_t &i, const hlf_t &
) { 
  return i-1; 
}
//listing06
template<class n_t>
inline rng_t operator^(
  const rng_t &r, const n_t &n
) { 
  return rng_t(
    (r - n).first(), 
    (r + n).last()
  ); 
} 
//listing07
template<int d> 
inline idx_2d_t pi(const rng_t &i, const rng_t &j);

template<>
inline idx_2d_t pi<0>(
  const rng_t &i, const rng_t &j
) {
  return idx_2d_t({i,j});
};

template<>
inline idx_2d_t pi<1>(
  const rng_t &j, const rng_t &i
) {
  return idx_2d_t({i,j});
}; 

#include "solvers.hpp"
#include "cyclic.hpp"
#include "donorcell.hpp"
#include "mpdata.hpp"
