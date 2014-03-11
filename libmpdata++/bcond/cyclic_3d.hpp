/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/bcond/cyclic_common.hpp>
#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template<int d, typename real_t>
    class cyclic_left_3d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void bcinit(const arr_t &a, const rng_t &j, const rng_t &k) {}

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
	a(pi<d>(this->left_halo_sclr, j, k)) = a(pi<d>(this->rght_intr_sclr, j, k)); 
      }

      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        a(pi<d>(this->left_halo_vctr, j, k)) = a(pi<d>(this->rght_intr_vctr, j, k));
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }
    };

    template<int d, typename real_t>
    class cyclic_rght_3d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void bcinit(const arr_t &a, const rng_t &j, const rng_t &k) {}

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
	a(pi<d>(this->rght_halo_sclr, j, k)) = a(pi<d>(this->left_intr_sclr, j, k));
      }
      
      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        a(pi<d>(this->rght_halo_vctr, j, k)) = a(pi<d>(this->left_intr_vctr, j, k));
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
