/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <advoocat/bcond/cyclic_common.hpp>

namespace advoocat
{
  namespace bcond
  {
    template <typename real_t>
    class cyclic_left_1d : public cyclic_left_common<real_t>
    {
      using parent_t = cyclic_left_common<real_t>;
      using arr_t = blitz::Array<real_t, 1>;

      public:

      // ctor
      cyclic_left_1d(const rng_t &i, int halo) : // TODO: inherit ctor
	parent_t(i, halo)
      { }

      // method invoked by the solver
      void fill_halos(const arr_t &a)
      {
	a(this->left_halo) = a(this->rght_edge);     
      }
    };

    template <typename real_t>
    class cyclic_rght_1d : public cyclic_rght_common<real_t>
    {
      using parent_t = cyclic_rght_common<real_t>;
      using arr_t = blitz::Array<real_t, 1>;

      public:

      // ctor
      cyclic_rght_1d(const rng_t &i, int halo) : // TODO: inherit ctor
	parent_t(i, halo)
      { }

      // method invoked by the solver
      void fill_halos(const arr_t &a)
      {
	a(this->rght_halo) = a(this->left_edge);     
      }
    };
  }; // namespace bcond
}; // namespace advoocat
