// 1D cyclic boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int dim>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     dim,
      typename std::enable_if<
        knd == cyclic &&
        dir == left   &&
        n_dims == 1
      >::type
    > : public detail::bcond_common<real_t, halo, n_dims>
    {
      using parent_t = detail::bcond_common<real_t, halo, n_dims>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(arr_t &a, const bool deriv = false)
      {
        a(this->left_halo_sclr) = a(this->rght_intr_sclr);
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const bool ad = false)
      {
        av[0](this->left_halo_vctr) = av[0](this->rght_intr_vctr);
      }

      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const bool ad = false)
      {
        fill_halos_vctr_alng(av, ad);
      }

      void copy_edge_sclr_to_halo1_cyclic(arr_t &a)
      {
        assert(halo>=1);

        a(this->left_halo_sclr.last()) = a(this->rght_edge_sclr);
      }

      void avg_edge_and_halo1_sclr_cyclic(arr_t &a)
      {
        assert(halo>=1);

        a(this->left_edge_sclr) = ( a(this->left_edge_sclr) + a(this->left_halo_sclr.last()) ) / real_t(2); 
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int dim>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     dim,
      typename std::enable_if<
        knd == cyclic &&
        dir == rght   &&
        n_dims == 1
      >::type
    > : public detail::bcond_common<real_t, halo, n_dims>
    {
      using parent_t = detail::bcond_common<real_t, halo, n_dims>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(arr_t &a, const bool deriv = false)
      {
        a(this->rght_halo_sclr) = a(this->left_intr_sclr);
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const bool ad = false)
      {
        av[0](this->rght_halo_vctr) = av[0](this->left_intr_vctr);
      }

      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const bool ad = false)
      {
        fill_halos_vctr_alng(av, ad);
      }

      void copy_edge_sclr_to_halo1_cyclic(arr_t &a)
      {
        assert(halo>=1);

        a(this->rght_halo_sclr.first()) = a(this->left_edge_sclr);
      }

      void avg_edge_and_halo1_sclr_cyclic(arr_t &a)
      {
        assert(halo>=1);

        a(this->rght_edge_sclr) = ( a(this->rght_edge_sclr) + a(this->rght_halo_sclr.first()) ) / real_t(2); 
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
