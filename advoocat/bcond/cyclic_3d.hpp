/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <advoocat/bcond/cyclic_common.hpp>
#include <advoocat/idxperm.hpp>

namespace advoocat
{
  namespace bcond
  {
    template<int d, typename real_t>
    class cyclic_left_3d : public cyclic_left_common<real_t>
    {
      using parent_t = cyclic_left_common<real_t>;
      using arr_t = blitz::Array<real_t, 3>;

      public:

      // ctor
      cyclic_left_3d(const rng_t &i, int halo) :
	parent_t(i, halo)
      {} 

      // method invoked by the solver
      void fill_halos(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
	a(pi<d>(this->left_halo, j, k)) = a(pi<d>(this->rght_edge, j, k));     
      }
    };

    template<int d, typename real_t>
    class cyclic_rght_3d : public cyclic_rght_common<real_t>
    {
      using parent_t = cyclic_rght_common<real_t>;
      using arr_t = blitz::Array<real_t, 3>;

      public:

      // ctor
      cyclic_rght_3d(const rng_t &i, int halo) :
	parent_t(i, halo)
      {} 

      // method invoked by the solver
      void fill_halos(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
	a(pi<d>(this->rght_halo, j, k)) = a(pi<d>(this->left_edge, j, k));     
      }
    };
  }; // namespace bcond
}; // namespace advoocat
