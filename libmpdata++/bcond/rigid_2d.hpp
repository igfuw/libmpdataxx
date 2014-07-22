/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template<int d, typename real_t>
    class rigid_left_2d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:
      
      void fill_halos_sclr(const arr_t &a, const rng_t &j, const bool deriv = false)
      {
        using namespace idxperm;
        // zero flux condition
        for (int i = this->left_halo_sclr.first(), n = this->halo; i <= this->left_halo_sclr.last(); ++i, --n)
        {
          a(pi<d>(i, j)) = a(pi<d>(this->left_edge_sclr + n, j));
        }
      }
      
      void fill_halos_pres(const arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        a(pi<d>(this->left_halo_sclr.last(), j)) = 2 * a(pi<d>(this->left_edge_sclr,     j))
                                                     - a(pi<d>(this->left_edge_sclr + 1, j));
      }
      
      void set_edge_pres(const arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        a(pi<d>(this->left_edge_sclr, j)) = 0;
      }
      
      void set_edge_pres(const arr_t &a, const arr_t &b, const rng_t &j)
      {
        using namespace idxperm;
        a(pi<d>(this->left_edge_sclr, j)) = -b(pi<d>(this->left_edge_sclr, j));
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
        // zero velocity condition
        for (int i = this->left_halo_vctr.first(), n = this->halo; i <= this->left_halo_vctr.last(); ++i, --n)
        {
	  av[d](pi<d>(i, j)) = -av[d](pi<d>(this->left_edge_sclr + n - h, j));
        }
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        // note intentional sclr
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
          a(pi<d>(i, j)) = 0; 
      }
    };

    template<int d, typename real_t>
    class rigid_rght_2d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor
      
      public:
      
      void fill_halos_sclr(const arr_t &a, const rng_t &j, const bool deriv = false)
      {
        // zero flux condition
        using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(), n = 1; i <= this->rght_halo_sclr.last(); ++i, ++n)
        {
          a(pi<d>(i, j)) = a(pi<d>(this->rght_edge_sclr - n, j)); // zero gradient for scalar gradient
        }
      }
      
      void fill_halos_pres(const arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        a(pi<d>(this->rght_halo_sclr.first(), j)) = 2 * a(pi<d>(this->rght_edge_sclr,     j))
                                                      - a(pi<d>(this->rght_edge_sclr - 1, j));
      }
      
      void set_edge_pres(const arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_edge_sclr, j)) = 0;
      }
      
      void set_edge_pres(const arr_t &a, const arr_t &b, const rng_t &j)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_edge_sclr, j)) = -b(pi<d>(this->rght_edge_sclr, j));
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
        // zero velocity condition
        for (int i = this->rght_halo_vctr.first(), n = 1; i <= this->rght_halo_vctr.last(); ++i, ++n)
        {
	  av[d](pi<d>(i, j)) = -av[d](pi<d>(this->rght_edge_sclr - n + h, j));
        }
      }
      
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        // note intentional sclr
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
          a(pi<d>(i, j)) = 0; 
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
