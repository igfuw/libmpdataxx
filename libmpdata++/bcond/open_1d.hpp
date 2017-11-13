// 1D open boundary conditions for libmpdata++
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
        knd == open &&
        dir == left &&
        n_dims == 1
      >::type
    > : public detail::bcond_common<real_t, halo>
    {
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor
      
      public:

      void fill_halos_sclr(arr_t &a, const bool deriv = false)
      {
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
        {
	  if (deriv)
            a(rng_t(i, i)) = 0;
          else
            a(rng_t(i, i)) = a(this->left_edge_sclr);
        }
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const bool ad = false)
      {
        for (int i = this->left_halo_vctr.first(); i <= this->left_halo_vctr.last() - (ad ? 1 : 0); ++i)
	  av[0](rng_t(i, i)) = av[0](this->left_intr_vctr.first());
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int dim>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     dim,
      typename std::enable_if<
        knd == open &&
        dir == rght &&
        n_dims == 1
      >::type
    > : public detail::bcond_common<real_t, halo>
    {
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor
      
      public:

      void fill_halos_sclr(arr_t &a, const bool deriv = false)
      {
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
        {
	  if (deriv) 
            a(rng_t(i, i)) = 0;
          else
            a(rng_t(i, i)) = a(this->rght_edge_sclr);
        }
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const bool ad = false)
      {
        for (int i = this->rght_halo_vctr.first() + (ad ? 1 : 0); i <= this->rght_halo_vctr.last(); ++i)
	  av[0](rng_t(i, i)) = av[0](this->rght_intr_vctr.first());
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
