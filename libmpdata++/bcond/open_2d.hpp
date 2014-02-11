/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/bcond/open_common.hpp>
#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template<int d, typename real_t>
    class open_left_2d : public open_left_common<real_t>
    {
      using parent_t = open_left_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
	a(pi<d>(this->left_halo_sclr, j)) = 0;
      }

      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        assert(parent_t::halo > 1 && "there is no vector halo for halo=1");
        a(pi<d>(this->left_halo_vctr, j)) = a(pi<d>(this->left_edge_vctr, j));
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        // note intentional sclr
        a(pi<d>(this->left_halo_sclr, j)) = a(pi<d>(this->left_edge_sclr, j));
      }
    };

    template<int d, typename real_t>
    class open_rght_2d : public open_rght_common<real_t>
    {
      using parent_t = open_rght_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
	a(pi<d>(this->rght_halo_sclr, j)) = 0;
      }

      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        assert(parent_t::halo > 1 && "there is no vector halo for halo=1");
        a(pi<d>(this->rght_halo_vctr, j)) = a(pi<d>(this->rght_edge_vctr, j));
      }
      
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        // note intentional sclr
        a(pi<d>(this->rght_halo_sclr, j)) = a(pi<d>(this->rght_edge_sclr, j));
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
