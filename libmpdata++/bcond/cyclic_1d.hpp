/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/bcond/cyclic_common.hpp>

namespace libmpdataxx
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

      // methods invoked by the solver
      void fill_halos_sclr(const arr_t &a)
      {
	a(this->left_halo_sclr) = a(this->rght_edge_sclr);
      }

      void fill_halos_vctr(const arr_t &a)
      {
        assert(parent_t::halo > 1 && "there is no vector halo for halo=1"); 
	a(this->left_halo_vctr) = a(this->rght_edge_vctr);
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

      // methods invoked by the solver
      void fill_halos_sclr(const arr_t &a)
      {
	a(this->rght_halo_sclr) = a(this->left_edge_sclr);
      }

      void fill_halos_vctr(const arr_t &a)
      {
        assert(parent_t::halo > 1 && "there is no vector halo for halo=1"); 
	a(this->rght_halo_vctr) = a(this->left_edge_vctr);
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
