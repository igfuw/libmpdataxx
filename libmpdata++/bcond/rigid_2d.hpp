// 2D rigid boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == rigid &&
        dir == left &&
        n_dims == 2
      >::type
    > : public detail::bcond_common<real_t, halo>
    { 
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:
      
      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
        using namespace idxperm;
        // zero flux condition
        for (int i = this->left_halo_sclr.first(), n = halo; i <= this->left_halo_sclr.last(); ++i, --n)
        {
          a(pi<d>(i, j)) = a(pi<d>(this->left_edge_sclr + n, j));
        }
      }
      
      void fill_halos_pres(arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        a(pi<d>(this->left_halo_sclr.last(), j)) = 2 * a(pi<d>(this->left_edge_sclr,     j))
                                                     - a(pi<d>(this->left_edge_sclr + 1, j));
      }

      void save_edge_vel(const arr_t &, const rng_t &) {}
     
      void set_edge_pres(arr_t &a, const rng_t &j, int)
      {
        using namespace idxperm;
        a(pi<d>(this->left_edge_sclr, j)) = 0;
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
	using namespace idxperm;
        // zero velocity condition
        for (int i = this->left_halo_vctr.first(), n = halo; i <= this->left_halo_vctr.last(); ++i, --n)
        {
	  av[d](pi<d>(i, j)) = -av[d](pi<d>(this->left_edge_sclr + n - h, j));
        }
      }

      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)
      {
        // note intentional sclr
        fill_halos_sclr(a, j);
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == rigid &&
        dir == rght &&
        n_dims == 2
      >::type
    > : public detail::bcond_common<real_t, halo>
    {
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor
      
      public:
      
      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
        // zero flux condition
        using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(), n = 1; i <= this->rght_halo_sclr.last(); ++i, ++n)
        {
          a(pi<d>(i, j)) = a(pi<d>(this->rght_edge_sclr - n, j)); // zero gradient for scalar gradient
        }
      }
      
      void fill_halos_pres(arr_t &a, const rng_t &j)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        a(pi<d>(this->rght_halo_sclr.first(), j)) = 2 * a(pi<d>(this->rght_edge_sclr,     j))
                                                      - a(pi<d>(this->rght_edge_sclr - 1, j));
      }
      
      void save_edge_vel(const arr_t &, const rng_t &) {}
      
      void set_edge_pres(arr_t &a, const rng_t &j, int)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_edge_sclr, j)) = 0;
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
	using namespace idxperm;
        // zero velocity condition
        for (int i = this->rght_halo_vctr.first(), n = 1; i <= this->rght_halo_vctr.last(); ++i, ++n)
        {
	  av[d](pi<d>(i, j)) = -av[d](pi<d>(this->rght_edge_sclr - n + h, j));
        }
      }
      
      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)
      {
        // note intentional sclr
        fill_halos_sclr(a, j);
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
