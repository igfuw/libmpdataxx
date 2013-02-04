/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "cyclic_common.hpp"
#include "../idxperm.hpp"

namespace advoocat
{
  namespace bcond
  {
// TODO: indent
template<int d, typename real_t = float>
class cyclic_3d : cyclic_common<3, real_t>
{
  using parent_t = cyclic_common<3, real_t>;
  using arr_3d_t = typename parent_t::arr_t;

  public:

  // ctor
  cyclic_3d(const rng_t &i, int halo) :
    parent_t(i, halo)
  {} 

  // method invoked by the solver
  void fill_halos(const arr_3d_t &a, const rng_t &j, const rng_t &k)
  {
    using namespace idxperm;
    a(pi<d>(this->left_halo, j, k)) = a(pi<d>(this->rght_edge, j, k));     
    a(pi<d>(this->rght_halo, j, k)) = a(pi<d>(this->left_edge, j, k));     
  }
};
  }; // namespace bcond
}; // namespace advoocat
