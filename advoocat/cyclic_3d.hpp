/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "cyclic_common.hpp"

template<int d, typename real_t = float>
class cyclic_3d : cyclic_common<real_t>
{
  public:

  // ctor
  cyclic_3d(const rng_t &i, int halo) :
    cyclic_common<real_t>(i, halo)
  {} 

  // method invoked by the solver
  void fill_halos(const arr_3d_t<real_t> &a, const rng_t &j, const rng_t &k)
  {
    a(pi<d>(this->left_halo, j, k)) = a(pi<d>(this->rght_edge, j, k));     
    a(pi<d>(this->rght_halo, j, k)) = a(pi<d>(this->left_edge, j, k));     
  }
};

