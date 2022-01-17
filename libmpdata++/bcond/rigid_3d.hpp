// 3D rigid boundary conditions for libmpdata++
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
        n_dims == 3
      >::type
    > : public detail::bcond_common<real_t, halo, n_dims>
    {
      using parent_t = detail::bcond_common<real_t, halo, n_dims>;
      using arr_t = blitz::Array<real_t, n_dims>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        using namespace idxperm;
        // zero flux condition
        for (int i = this->left_halo_sclr.first(), n = halo; i <= this->left_halo_sclr.last(); ++i, --n)
        {
          a(pi<d>(i, j, k)) = a(pi<d>(this->left_edge_sclr + n, j, k));
        }
      }

      void fill_halos_pres(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        for (int i = this->left_halo_sclr.first(), n = halo; i <= this->left_halo_sclr.last(); ++i, --n)
        {
          a(pi<d>(i, j, k)) = 2 * a(pi<d>(this->left_edge_sclr,     j, k))
                                - a(pi<d>(this->left_edge_sclr + n, j, k));
        }
      }

      void save_edge_vel(const arr_t &, const rng_t &, const rng_t &) {}

      void set_edge_pres(arr_t &a, const rng_t &j, const rng_t &k, int)
      {
        using namespace idxperm;
        a(pi<d>(this->left_edge_sclr, j, k)) = 0;
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
        using namespace idxperm;
        // zero velocity condition
        for (int i = this->left_halo_vctr.first(), n = halo; i <= this->left_halo_vctr.last() - (ad ? 1 : 0); ++i, --n)
        {
          av[d](pi<d>(i, j, k)) = -av[d](pi<d>(this->left_edge_sclr + n - h, j, k));
        }
      }

      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j, const rng_t &k)
      {
        // note intentional sclr
        fill_halos_sclr(a, j, k);
      }

      void fill_halos_flux(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // zero flux condition
        av[d](pi<d>(this->left_halo_vctr.last(), j, k)) = -av[d](pi<d>(this->left_edge_sclr + h, j, k));
      }

      void fill_halos_sgs_div(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        a(pi<d>(this->left_edge_sclr - h, j, k)) = 2 * a(pi<d>(this->left_edge_sclr + h, j, k))
                                                   -   a(pi<d>(this->left_edge_sclr + 1 + h, j, k));
      }

      void fill_halos_sgs_div_stgr(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sgs_div(a, j, k);
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const rng_t &j, const rng_t &k, const int offset = 0)
      {
        // fill halos for a staggered field so that it has zero value on the edge
        // that is 0.5 * (a(edge-h) + a(edge+h)) = 0
        using namespace idxperm;
        const auto &a = av[offset + d];
        a(pi<d>(this->left_edge_sclr - h, j, k)) = - a(pi<d>(this->left_edge_sclr + h, j, k));
      }

      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const rng_t &j, const rng_t &k, const real_t di)
      {
        using namespace idxperm;
        const auto &a = av[d];
        a(pi<d>(this->left_edge_sclr - h, j, k)) = 2 * a(pi<d>(this->left_edge_sclr + h, j, k))
                                                   -   a(pi<d>(this->left_edge_sclr + 1 + h, j, k));
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == rigid &&
        dir == rght &&
        n_dims == 3
      >::type
    > : public detail::bcond_common<real_t, halo, n_dims>
    {
      using parent_t = detail::bcond_common<real_t, halo, n_dims>;
      using arr_t = blitz::Array<real_t, n_dims>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        // zero flux condition
        using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(), n = 1; i <= this->rght_halo_sclr.last(); ++i, ++n)
        {
          a(pi<d>(i, j, k)) = a(pi<d>(this->rght_edge_sclr - n, j, k)); // zero gradient for scalar gradient
        }
      }


      void save_edge_vel(const arr_t &, const rng_t &, const rng_t &) {}

      void fill_halos_pres(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        for (int i = this->rght_halo_sclr.first(), n = 1; i <= this->rght_halo_sclr.last(); ++i, ++n)
        {
          a(pi<d>(i, j, k)) = 2 * a(pi<d>(this->rght_edge_sclr,     j, k))
                                - a(pi<d>(this->rght_edge_sclr - n, j, k));
        }
      }

      void set_edge_pres(arr_t &a, const rng_t &j, const rng_t &k, int)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_edge_sclr, j, k)) = 0;
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
        using namespace idxperm;
        // zero velocity condition
        for (int i = this->rght_halo_vctr.first() + (ad ? 1 : 0), n = 1; i <= this->rght_halo_vctr.last(); ++i, ++n)
        {
          av[d](pi<d>(i, j, k)) = -av[d](pi<d>(this->rght_edge_sclr - n + h, j, k));
        }
      }

      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j, const rng_t &k)
      {
        // note intentional sclr
        fill_halos_sclr(a, j, k);
      }

      void fill_halos_flux(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // zero flux condition
        av[d](pi<d>(this->rght_halo_vctr.first(), j, k)) = -av[d](pi<d>(this->rght_edge_sclr - h, j, k));
      }

      void fill_halos_sgs_div(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_edge_sclr + h, j, k)) = 2 * a(pi<d>(this->rght_edge_sclr - h, j, k))
                                                   -   a(pi<d>(this->rght_edge_sclr - 1 - h, j, k));
      }

      void fill_halos_sgs_div_stgr(arr_t &a, const rng_t &j, const rng_t &k)
      {
        fill_halos_sgs_div(a, j, k);
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const rng_t &j, const rng_t &k, const int offset = 0)
      {
        // fill halos for a staggered field so that it has zero value on tke edge
        // that is 0.5 * (a(edge-h) + a(edge+h)) = 0
        using namespace idxperm;
        const auto &a = av[offset + d];
        a(pi<d>(this->rght_edge_sclr + h, j, k)) = - a(pi<d>(this->rght_edge_sclr - h, j, k));
      }

      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const rng_t &j, const rng_t &k, const real_t di)
      {
        using namespace idxperm;
        const auto &a = av[d];
        a(pi<d>(this->rght_edge_sclr + h, j, k)) = 2 * a(pi<d>(this->rght_edge_sclr - h, j, k))
                                                   -   a(pi<d>(this->rght_edge_sclr - 1 - h, j, k));
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
