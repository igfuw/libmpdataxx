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
    class open_left_3d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      real_t init_sclr;

      public:

      // method invoked by the solver
      void bcinit(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        assert(min(a(pi<d>(this->left_edge_sclr, j, k))) == max(a(pi<d>(this->left_edge_sclr, j, k)))
               && "variable initial signal on open boundary");
	init_sclr = min(a(pi<d>(this->left_edge_sclr, j, k)));
      }

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
	  a(pi<d>(rng_t(i, i), j, k)) = init_sclr;
      }

      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        for (int i = this->left_halo_vctr.first(); i <= this->left_halo_vctr.last(); ++i)
          a(pi<d>(rng_t(i, i), j, k)) = a(pi<d>(rng_t(this->left_intr_vctr.first(),
                                                   this->left_intr_vctr.first()),
                                             j, k));
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
          a(pi<d>(rng_t(i, i), j, k)) = a(pi<d>(this->left_edge_sclr, j, k));
      }
    };

    template<int d, typename real_t>
    class open_rght_3d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor
      
      real_t init_sclr;

      public:

      // method invoked by the solver
      void bcinit(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        assert(min(a(pi<d>(this->rght_edge_sclr, j, k))) == max(a(pi<d>(this->rght_edge_sclr, j, k)))
               && "variable initial signal on open boundary");
	init_sclr = min(a(pi<d>(this->rght_edge_sclr, j, k)));
      }

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
	  a(pi<d>(rng_t(i, i), j, k)) = init_sclr;
      }

      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        for (int i = this->rght_halo_vctr.first(); i <= this->rght_halo_vctr.last(); ++i)
          a(pi<d>(rng_t(i, i), j, k)) = a(pi<d>(rng_t(this->rght_intr_vctr.last(),
                                                   this->rght_intr_vctr.last()),
                                             j, k));
      }
      
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
          a(pi<d>(rng_t(i, i), j, k)) = a(pi<d>(this->rght_edge_sclr, j, k));
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
