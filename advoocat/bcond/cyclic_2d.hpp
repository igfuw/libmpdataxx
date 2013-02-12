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
    template<int d, typename real_t = float>
    class cyclic_2d : public cyclic_common<2, real_t>
    {
      using parent_t = cyclic_common<2, real_t>;
      using arr_2d_t = typename parent_t::arr_t;

      public:

      // ctor
      cyclic_2d(const rng_t &i, int halo) :
	parent_t(i, halo)
      {} 

      // method invoked by the solver
      void fill_halos(const arr_2d_t &a, const rng_t &j)
      {
	using namespace idxperm;
	a(pi<d>(this->left_halo, j)) = a(pi<d>(this->rght_edge, j));     
	a(pi<d>(this->rght_halo, j)) = a(pi<d>(this->left_edge, j));     
      }
    };
  }; // namespace bcond
}; // namespace advoocat
