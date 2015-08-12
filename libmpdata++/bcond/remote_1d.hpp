// 1D MPI ``remote'' boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/remote_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int dim>    
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     dim,
      typename std::enable_if<
        knd == remote && 
        dir == left   && 
        n_dims == 1
      >::type
    > : public detail::remote_common<real_t, halo, dir>
    {
      using parent_t = detail::remote_common<real_t, halo, dir>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : -1;

      public:

      void fill_halos_sclr(const arr_t &a, const bool deriv = false)
      {
        this->xchng(a, this->left_intr_sclr + off, this->left_halo_sclr);
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av)
      {
        this->xchng(av[0], this->left_intr_vctr + off, this->left_halo_vctr);
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int dim>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     dim,
      typename std::enable_if<
        knd == remote &&
        dir == rght   &&
        n_dims == 1
      >::type
    > : public detail::remote_common<real_t, halo, dir>
    {
      using parent_t = detail::remote_common<real_t, halo, dir>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : 1;

      public:

      void fill_halos_sclr(const arr_t &a, const bool deriv = false)
      {
        this->xchng(a, this->rght_intr_sclr + off, this->rght_halo_sclr);
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av)
      {
        this->xchng(av[0], this->rght_intr_vctr + off, this->rght_halo_vctr);
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
