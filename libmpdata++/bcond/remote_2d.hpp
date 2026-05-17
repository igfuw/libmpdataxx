// 2D MPI ``remote'' boundary conditions for libmpdata++
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
        n_dims == 2
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {
      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : -1;

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
        using namespace idxperm;
        this->xchng(a, pi<d>(this->left_intr_sclr + off, j), pi<d>(this->left_halo_sclr, j));
      }

      void fill_halos_pres(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }

      void save_edge_vel(const arr_t &, const rng_t &) {}

      void save_edge_val(const arr_t &, const arr_t &, const rng_t &) {}

      void set_edge_pres(arr_t &, const rng_t &, int) {}

      // we require that given process calculates its internal vectors + the nearest vector to the left, TODO: how to enforce this?
      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
        using namespace idxperm;
        if(!this->is_cyclic)
        {
          if(halo == 1)
            // send vectors to the left of the domain
            this->send(av[0], pi<d>(this->left_intr_vctr + off, j));
          else
            // receive the halo without the rightmost column, which was caluclated by this process
            this->xchng(av[0], pi<d>(this->left_intr_vctr + off, j), pi<d>((this->left_halo_vctr^h)^(-1), j));
        }
        else
          this->xchng(av[0], pi<d>(this->left_intr_vctr + off, j), pi<d>(this->left_halo_vctr, j));
      }

      void fill_halos_sgs_div(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }

      void fill_halos_sgs_div_stgr(arr_t &a, const rng_t &j)
      {
        fill_halos_sgs_div(a, j); // NOTE: probably should be replaced with fill_halos_vctr_alng, but not an issue right now because there is no remote bcond in the vertical and vip_div is staggered only in the vertical (?)
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const rng_t &j, const int offset = 0)
      {
        using namespace idxperm;
        // the same logic as fill_halos_vctr_alng but have to consider offset ... TODO: find a way to reuse !
        if(!this->is_cyclic)
        {
          if(halo == 1)
            // send vectors to the left of the domain
            this->send(av[0 + offset], pi<d>(this->left_intr_vctr + off, j));
          else
            // receive the halo without the rightmost column, which was caluclated by this process
            this->xchng(av[0 + offset], pi<d>(this->left_intr_vctr + off, j), pi<d>((this->left_halo_vctr^h)^(-1), j));
        }
        else
          this->xchng(av[0 + offset], pi<d>(this->left_intr_vctr + off, j), pi<d>(this->left_halo_vctr, j));
      }

      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const rng_t &j, const real_t)
      {
        fill_halos_vctr_alng(av, j);
      }

      // TODO: move to common? (same in cyclic!)
      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }

      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
        fill_halos_vctr_alng(av, j, ad);
      }

      void fill_halos_vctr_nrml_cyclic(arr_t &a, const rng_t &j)
      {
        fill_halos_vctr_nrml(a, j);
      }

      void copy_edge_sclr_to_halo1_cyclic(arr_t &a, const rng_t &j)
      {
        if(!this->is_cyclic)
          return;

        using namespace idxperm;
        assert(halo>=1);
        this->xchng(a, pi<d>(this->left_edge_sclr, j), pi<d>(this->left_halo_sclr.last(), j));
      }

      void avg_edge_and_halo1_sclr_cyclic(arr_t &a, const rng_t &j)
      {
        if(!this->is_cyclic)
          return;

        using namespace idxperm;
        assert(halo>=1);
        a(pi<d>(this->left_edge_sclr, j)) = ( a(pi<d>(this->left_edge_sclr, j)) + a(pi<d>(this->left_halo_sclr.last(), j)) ) / real_t(2);
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == remote &&
        dir == rght   &&
        n_dims == 2
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {
      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : 1;

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
        using namespace idxperm;
        this->xchng(a, pi<d>(this->rght_intr_sclr + off, j), pi<d>(this->rght_halo_sclr, j));
      }

      void fill_halos_pres(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }

      void save_edge_vel(const arr_t &, const rng_t &) {}

      void save_edge_val(const arr_t &, const arr_t &, const rng_t &) {}

      void set_edge_pres(arr_t &, const rng_t &, int) {}

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
        using namespace idxperm;
        if(!this->is_cyclic)
        {
          if(halo == 1)
            //receive the halo
            this->recv(av[0], pi<d>(this->rght_halo_vctr, j));
          else
            // don't send the first column to the right of the domain, it will be calculated and sent here by the process to the right
            this->xchng(av[0], pi<d>(((this->rght_intr_vctr + off)^h)^(-1), j), pi<d>(this->rght_halo_vctr, j));
        }
        else
          this->xchng(av[0], pi<d>(this->rght_intr_vctr + off, j), pi<d>(this->rght_halo_vctr, j));
      }

      void fill_halos_sgs_div(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }

      void fill_halos_sgs_div_stgr(arr_t &a, const rng_t &j)
      {
        fill_halos_sgs_div(a, j);
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const rng_t &j, const int offset = 0)
      {
        using namespace idxperm;
        // the same logic as fill_halos_vctr_alng but have to consider offset ... TODO: find a way to reuse !
        if(!this->is_cyclic)
        {
          if(halo == 1)
            //receive the halo
            this->recv(av[0 + offset], pi<d>(this->rght_halo_vctr, j));
          else
            // don't send the first column to the right of the domain, it will be calculated and sent here by the process to the right
            this->xchng(av[0 + offset], pi<d>(((this->rght_intr_vctr + off)^h)^(-1), j), pi<d>(this->rght_halo_vctr, j));
        }
        else
          this->xchng(av[0 + offset], pi<d>(this->rght_intr_vctr + off, j), pi<d>(this->rght_halo_vctr, j));
      }

      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const rng_t &j, const real_t)
      {
        fill_halos_vctr_alng(av, j);
      }

      // TODO: move to common? (same in cyclic!)
      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }

      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
        fill_halos_vctr_alng(av, j, ad);
      }

      void fill_halos_vctr_nrml_cyclic(arr_t &a, const rng_t &j)
      {
        fill_halos_vctr_nrml(a, j);
      }

      void copy_edge_sclr_to_halo1_cyclic(arr_t &a, const rng_t &j)
      {
        if(!this->is_cyclic)
          return;

        using namespace idxperm;
        assert(halo>=1);
        this->xchng(a, pi<d>(this->rght_edge_sclr, j), pi<d>(this->rght_halo_sclr.first(), j));
      }

      void avg_edge_and_halo1_sclr_cyclic(arr_t &a, const rng_t &j)
      {
        if(!this->is_cyclic)
          return;

        using namespace idxperm;
        assert(halo>=1);
        a(pi<d>(this->rght_edge_sclr, j)) = ( a(pi<d>(this->rght_edge_sclr, j)) + a(pi<d>(this->rght_halo_sclr.first(), j)) ) / real_t(2);
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
