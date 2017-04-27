// 3D MPI ``remote'' boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/remote_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>    
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == remote && 
        dir == left   && 
        n_dims == 3
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {

      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : -1;

      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        using namespace idxperm;
        this->xchng(a, pi<d>(this->left_intr_sclr + off, j, k), pi<d>(this->left_halo_sclr, j, k));
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        this->xchng(av[0], pi<d>(this->left_intr_vctr + off, j, k), pi<d>(this->left_halo_vctr, j, k));
      }

      // TODO: move to common? (same in cyclic!)
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)                 
      {                                                                         
        fill_halos_sclr(a, j, k);                                                  
      }  

    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == remote &&
        dir == rght   &&
        n_dims == 3
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {
      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : 1;

      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        using namespace idxperm;
        this->xchng(a, pi<d>(this->rght_intr_sclr + off, j, k), pi<d>(this->rght_halo_sclr, j, k));
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        this->xchng(av[0], pi<d>(this->rght_intr_vctr + off, j, k), pi<d>(this->rght_halo_vctr, j, k));
      }

      // TODO: move to common? (same in cyclic!)
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)                 
      {                                                                         
        fill_halos_sclr(a, j, k);                                                  
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
