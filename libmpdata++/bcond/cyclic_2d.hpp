/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/bcond/cyclic_common.hpp>
#include <libmpdata++/idxperm.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template<int d, typename real_t>
    class cyclic_left_2d : public cyclic_left_common<real_t>
    {
      using parent_t = cyclic_left_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;

      public:

      // ctor
      cyclic_left_2d(const rng_t &i, int halo) : // TODO: inherit ctor
	parent_t(i, halo)
      {} 

      // method invoked by the solver
      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
	a(pi<d>(this->left_halo_sclr, j)) = a(pi<d>(this->rght_edge_sclr, j));
      }

      void fill_halos_vctr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        assert(parent_t::halo > 1 && "there is no vector halo for halo=1");
        a(pi<d>(this->left_halo_vctr, j)) = a(pi<d>(this->rght_edge_vctr, j));
      }
    };

    template<int d, typename real_t>
    class cyclic_rght_2d : public cyclic_rght_common<real_t>
    {
      using parent_t = cyclic_rght_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;

      public:

      // ctor
      cyclic_rght_2d(const rng_t &i, int halo) : // TODO: inherit ctor
	parent_t(i, halo)
      {} 

      // method invoked by the solver
      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
	a(pi<d>(this->rght_halo_sclr, j)) = a(pi<d>(this->left_edge_sclr, j));
      }

      void fill_halos_vctr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        assert(parent_t::halo > 1 && "there is no vector halo for halo=1");
        a(pi<d>(this->rght_halo_vctr, j)) = a(pi<d>(this->left_edge_vctr, j));
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
