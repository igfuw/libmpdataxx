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

struct cyclic_1d
{
  // member fields
  idx_1d_t left_halo, rght_halo;
  idx_1d_t left_edge, rght_edge;;

  // ctor
  cyclic_1d(
    const rng_t &i, int hlo
  ) :
    left_halo( rng_t(i.first()-hlo,   i.first()-1    )),
    rght_edge( rng_t(i.last() -hlo+1, i.last()       )),
    rght_halo( rng_t(i.last() +1,     i.last() +hlo  )),
    left_edge( rng_t(i.first(),       i.first()+hlo-1))
  {} 

  // method invoked by the solver
  void fill_halos(const arr_1d_t &a)
  {
    a(left_halo) = a(rght_edge);     
    a(rght_halo) = a(left_edge);     
  }
};

template<int d>
struct cyclic_2d
{
  // member fields
  idx_2d_t left_halo, rght_halo;
  idx_2d_t left_edge, rght_edge;;

  // ctor
  cyclic_2d(
    const rng_t &i, const rng_t &j, int hlo
  ) :
    left_halo(pi<d>(
      rng_t(i.first()-hlo, i.first()-1), j//^hlo // TODO Range::all()
    )),
    rght_edge(pi<d>(
      rng_t(i.last()-hlo+1, i.last()  ), j//^hlo // TODO Range::all()
    )),
    rght_halo(pi<d>(
      rng_t(i.last()+1, i.last()+hlo  ), j//^hlo // TODO Range::all()
    )),
    left_edge(pi<d>(
      rng_t(i.first(), i.first()+hlo-1), j//^hlo // TODO Range::all()
    ))
  {} 

  // method invoked by the solver
  void fill_halos(const arr_2d_t &a)
  {
    a(left_halo) = a(rght_edge);     
    a(rght_halo) = a(left_edge);     
  }
};

#include "donorcell.hpp"
#include "mpdata.hpp"
