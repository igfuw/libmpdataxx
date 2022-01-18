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
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == cyclic &&
        dir == left &&
        n_dims == 3
      >::type
    > : public detail::bcond_common<real_t, halo, n_dims>
    {
      using parent_t = detail::bcond_common<real_t, halo, n_dims>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        using namespace idxperm;
        a(pi<d>(this->left_halo_sclr, j, k)) = a(pi<d>(this->rght_intr_sclr, j, k));
      }

      void fill_halos_pres(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }

      void save_edge_vel(const arr_t &, const rng_t &, const rng_t &) {}

      void set_edge_pres(arr_t &, const rng_t &, const rng_t &, int) {}

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
        using namespace idxperm;
        av[d](pi<d>(this->left_halo_vctr, j, k)) = av[d](pi<d>(this->rght_intr_vctr, j, k));
      }

      void fill_halos_sgs_div(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }

      void fill_halos_sgs_div_stgr(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        a(pi<d>(this->left_halo_vctr, j, k)) = a(pi<d>(this->rght_intr_vctr, j, k));
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const rng_t &j, const rng_t &k, const int offset = 0)
      {
        using namespace idxperm;
        // the same logic as fill_halos_vctr_alng but have to consider offset ... TODO: find a way to reuse !
        av[d + offset](pi<d>(this->left_halo_vctr, j, k)) = av[d + offset](pi<d>(this->rght_intr_vctr, j, k));
      }

      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const rng_t &j, const rng_t &k, const real_t)
      {
        fill_halos_vctr_alng(av, j, k);
      }

      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }

      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
        fill_halos_vctr_alng(av, j, k, ad);
      }

      void fill_halos_vctr_nrml_cyclic(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_vctr_nrml(a, j, k);
      }

      void copy_edge_sclr_to_halo1_cyclic(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        assert(halo>=1);

        a(pi<d>(this->left_halo_sclr.last(), j, k)) = a(pi<d>(this->rght_edge_sclr, j, k));
      }

      void avg_edge_and_halo1_sclr_cyclic(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        assert(halo>=1);

        a(pi<d>(this->left_edge_sclr, j, k)) = ( a(pi<d>(this->left_edge_sclr, j, k)) + a(pi<d>(this->left_halo_sclr.last(), j, k)) ) / real_t(2); 
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == cyclic &&
        dir == rght &&
        n_dims == 3
      >::type
    > : public detail::bcond_common<real_t, halo, n_dims>
    {
      using parent_t = detail::bcond_common<real_t, halo, n_dims>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_halo_sclr, j, k)) = a(pi<d>(this->left_intr_sclr, j, k));
      }

      void fill_halos_pres(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }

      void save_edge_vel(const arr_t &, const rng_t &, const rng_t &) {}

      void set_edge_pres(arr_t &, const rng_t &, const rng_t &, int) {}

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
        using namespace idxperm;
        av[d](pi<d>(this->rght_halo_vctr, j, k)) = av[d](pi<d>(this->left_intr_vctr, j, k));
      }

      void fill_halos_sgs_div(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }

      void fill_halos_sgs_div_stgr(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_halo_vctr, j, k)) = a(pi<d>(this->left_intr_vctr, j, k));
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const rng_t &j, const rng_t &k, const int offset = 0)
      {
        using namespace idxperm;
        // the same logic as fill_halos_vctr_alng but have to consider offset ... TODO: find a way to reuse !
        av[d + offset](pi<d>(this->rght_halo_vctr, j, k)) = av[d + offset](pi<d>(this->left_intr_vctr, j, k));
      }

      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const rng_t &j, const rng_t &k, const real_t)
      {
        fill_halos_vctr_alng(av, j, k);
      }

      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sclr(a, j, k);
      }

      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
        fill_halos_vctr_alng(av, j, k, ad);
      }

      void fill_halos_vctr_nrml_cyclic(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_vctr_nrml(a, j, k);
      }

      void copy_edge_sclr_to_halo1_cyclic(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        assert(halo>=1);

        a(pi<d>(this->rght_halo_sclr.first(), j, k)) = a(pi<d>(this->left_edge_sclr, j, k));
      }

      void avg_edge_and_halo1_sclr_cyclic(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        assert(halo>=1);

        a(pi<d>(this->rght_edge_sclr, j, k)) = ( a(pi<d>(this->rght_edge_sclr, j, k)) + a(pi<d>(this->rght_halo_sclr.first(), j, k)) ) / real_t(2); 
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
