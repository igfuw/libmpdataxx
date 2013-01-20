/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "cyclic_common.hpp"

template <typename real_t = float>
class cyclic_1d : cyclic_common<real_t>
{
  public:

  // ctor
  cyclic_1d(const rng_t &i, int halo) :
    cyclic_common<real_t>(i, halo)
  {}

  // method invoked by the solver
  void fill_halos(const arr_1d_t<real_t> &a)
  {
    a(this->left_halo) = a(this->rght_edge);     
    a(this->rght_halo) = a(this->left_edge);     
  }
};
