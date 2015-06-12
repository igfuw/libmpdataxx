// 3D cyclic boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == cyclic &&
        dir == left &&
        n_dims == 3
      >::type
    > : public detail::bcond_common<real_t>
    { 
      using parent_t = detail::bcond_common<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
	using namespace idxperm;
	a(pi<d>(this->left_halo_sclr, j, k)) = a(pi<d>(this->rght_intr_sclr, j, k)); 
      }

      void fill_halos_pres(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }
      
      void set_edge_pres(const arr_t &a, const rng_t &j, const rng_t &k) {}
      
      void set_edge_pres(const arr_t &a, const arr_t &b, const rng_t &j, const rng_t &k) {}

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        av[d](pi<d>(this->left_halo_vctr, j, k)) = av[d](pi<d>(this->rght_intr_vctr, j, k));
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }
    };

    template <typename real_t, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == cyclic &&
        dir == rght &&
        n_dims == 3
      >::type
    > : public detail::bcond_common<real_t>
    {
      using parent_t = detail::bcond_common<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
	using namespace idxperm;
	a(pi<d>(this->rght_halo_sclr, j, k)) = a(pi<d>(this->left_intr_sclr, j, k));
      }

      void fill_halos_pres(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }
      
      void set_edge_pres(const arr_t &a, const rng_t &j, const rng_t &k) {}

      void set_edge_pres(const arr_t &a, const arr_t &b, const rng_t &j, const rng_t &k) {}
      
      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        av[d](pi<d>(this->rght_halo_vctr, j, k)) = av[d](pi<d>(this->left_intr_vctr, j, k));
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
