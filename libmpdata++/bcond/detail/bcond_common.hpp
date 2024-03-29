// common code for all boundary conditions
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/formulae/arakawa_c.hpp>
#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    using namespace arakawa_c;

    enum bcond_e { null, cyclic, polar, open, rigid, remote, gndsky, custom };
    enum drctn_e { left, rght };

    template<
      typename real_t,
      int halo,
      bcond_e knd,
      drctn_e dir,
      int n_dims,
      int dim,
      class enableif = void
    >
    class bcond
    {};

    namespace detail
    {
      template <typename real_t, int halo, int n_dims>
      class bcond_common
      {
        protected:

        using arr_1d_t = blitz::Array<real_t, 1>;
        using arr_2d_t = blitz::Array<real_t, 2>;
        using arr_3d_t = blitz::Array<real_t, 3>;

        // needed by remote_3d and open_3d, stored here for convenience and potential use in other bconds
        const int thread_rank, thread_size;

        public:

        // 1D
        virtual void fill_halos_sclr(arr_1d_t &, const bool deriv = false)
        {
          assert(false && "bcond::fill_halos_sclr() called!");
        };

        virtual void fill_halos_vctr_alng(arrvec_t<arr_1d_t> &, const bool ad = false)
        {
          assert(false && "bcond::fill_halos_vctr() called!");
        };

        virtual void fill_halos_vctr_alng_cyclic(arrvec_t<blitz::Array<real_t, 1>> &, const bool ad = false)
        {};

        virtual void copy_edge_sclr_to_halo1_cyclic(arr_1d_t &)
        {};

        virtual void avg_edge_and_halo1_sclr_cyclic(arr_1d_t &)
        {};

        // 2D
        virtual void fill_halos_sclr(arr_2d_t &, const rng_t &, const bool deriv = false)
        {
          assert(false && "bcond::fill_halos_sclr() called!");
        };

        virtual void fill_halos_pres(arr_2d_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_pres() called!");
        };

        virtual void save_edge_vel(const arr_2d_t &, const rng_t &)
        {
          assert(false && "bcond::save_edge_vel() called!");
        };

        virtual void set_edge_pres(arr_2d_t &, const rng_t &, int)
        {
          assert(false && "bcond::set_edge() called!");
        };

        virtual void fill_halos_vctr_alng(arrvec_t<arr_2d_t> &, const rng_t &, const bool ad = false)
        {
          assert(false && "bcond::fill_halos_vctr_alng() called!");
        };

        virtual void fill_halos_sgs_div(arr_2d_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_sgs_div() called!");
        };

        virtual void fill_halos_sgs_div_stgr(arr_2d_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_sgs_div_stgr() called!");
        };

        virtual void fill_halos_sgs_vctr(arrvec_t<arr_2d_t> &, const arr_2d_t &, const rng_t &, const int offset = 0)
        {
          assert(false && "bcond::fill_halos_sgs_vctr() called!");
        };

        virtual void fill_halos_sgs_tnsr(arrvec_t<arr_2d_t> &, const arr_2d_t &, const arr_2d_t &, const rng_t &, const real_t)
        {
          assert(false && "bcond::fill_halos_sgs_tnsr called!");
        };

        virtual void fill_halos_vctr_nrml(arr_2d_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_vctr_nrml() called!");
        };

        virtual void fill_halos_vctr_alng_cyclic(arrvec_t<blitz::Array<real_t, 2>> &, const rng_t &, const bool ad = false)
        {};

        virtual void fill_halos_vctr_nrml_cyclic(blitz::Array<real_t, 2> &, const rng_t &)
        {};

        virtual void fill_halos_flux(arrvec_t<blitz::Array<real_t, 2>> &, const rng_t &)
        {};

        virtual void copy_edge_sclr_to_halo1_cyclic(arr_2d_t &, const rng_t &)
        {};

        virtual void avg_edge_and_halo1_sclr_cyclic(arr_2d_t &, const rng_t &)
        {};

        // 3D
        virtual void fill_halos_sclr(arr_3d_t &, const rng_t &, const rng_t &, const bool deriv = false)
        {
          assert(false && "bcond::fill_halos_sclr() called!");
        };

        virtual void fill_halos_pres(arr_3d_t &, const rng_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_pres() called!");
        };

        virtual void save_edge_vel(const arr_3d_t &, const rng_t &, const rng_t &)
        {
          assert(false && "bcond::save_edge_vel() called!");
        };

        virtual void set_edge_pres(arr_3d_t &, const rng_t &, const rng_t &, int)
        {
          assert(false && "bcond::set_edge() called!");
        };

        virtual void fill_halos_vctr_alng(arrvec_t<arr_3d_t> &, const rng_t &, const rng_t &, const bool ad = false)
        {
          assert(false && "bcond::fill_halos_vctr() called!");
        };

        virtual void fill_halos_sgs_div(arr_3d_t &, const rng_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_sgs_div() called!");
        };

        virtual void fill_halos_sgs_div_stgr(arr_3d_t &, const rng_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_sgs_div_stgr() called!");
        };

        virtual void fill_halos_sgs_vctr(arrvec_t<arr_3d_t> &,
                                            const arr_3d_t &,
                                            const rng_t &,
                                            const rng_t &,
                                            const int offset = 0)
        {
          assert(false && "bcond::fill_halos_sgs_vctr() called!");
        };

        virtual void fill_halos_sgs_tnsr(arrvec_t<arr_3d_t> &,
                                            const arr_3d_t &,
                                            const arr_3d_t &,
                                            const rng_t &,
                                            const rng_t &,
                                            const real_t)
        {
          assert(false && "bcond::fill_halos_sgs_tnsr called!");
        };

        virtual void fill_halos_vctr_nrml(arr_3d_t &, const rng_t &, const rng_t &)
        {
          assert(false && "bcond::fill_halos_vctr_nrml() called!");
        };

        virtual void fill_halos_vctr_alng_cyclic(arrvec_t<blitz::Array<real_t, 3>> &, const rng_t &, const rng_t &, const bool ad = false)
        {};

        virtual void fill_halos_vctr_nrml_cyclic(blitz::Array<real_t, 3> &, const rng_t &, const rng_t &)
        {};

        virtual void fill_halos_flux(arrvec_t<blitz::Array<real_t, 3>> &, const rng_t &, const rng_t &)
        {};

        virtual void copy_edge_sclr_to_halo1_cyclic(arr_3d_t &, const rng_t &, const rng_t &)
        {};

        virtual void avg_edge_and_halo1_sclr_cyclic(arr_3d_t &, const rng_t &, const rng_t &)
        {};

        const bool single_threaded;

        protected:
          // sclr
        int
          left_edge_sclr, rght_edge_sclr;
        rng_t
          left_halo_sclr, rght_halo_sclr,
          left_intr_sclr, rght_intr_sclr,
          // vctr
          left_halo_vctr, rght_halo_vctr,
          left_intr_vctr, rght_intr_vctr;

        public:

        // ctor
        bcond_common(
          const rng_t &i,
          const std::array<int, n_dims> &,
          bool single_threaded = false,
          const int thread_rank = -1, // -1 to indicate undefined
          const int thread_size = -1  // ditto
        ) :
          // sclr
          left_edge_sclr(
            i.first()
          ),
          rght_edge_sclr(
            i.last()
          ),
          left_halo_sclr(
            (i^halo).first(),
            (i^halo).first() + halo - 1
          ),
          rght_halo_sclr(
            (i^halo).last() - (halo - 1),
            (i^halo).last()
          ),
          left_intr_sclr(
            (i^(-1)).first(),
            (i^(-1)).first() + halo - 1
          ),
          rght_intr_sclr(
            (i^(-1)).last() - (halo - 1),
            (i^(-1)).last()
          ),
          // vctr
          left_halo_vctr(
            (i^h^(halo-1)).first(),
            (i^h^(halo-1)).first() + halo - 1
          ),
          rght_halo_vctr(
            (i^h^(halo-1)).last() - (halo - 1),
            (i^h^(halo-1)).last()
          ),
          left_intr_vctr(
            (i^h^(-1)).first(),
            (i^h^(-1)).first() + halo - 1
          ),
          rght_intr_vctr(
            (i^h^(-1)).last() - (halo - 1),
            (i^h^(-1)).last()
          ),
          single_threaded(single_threaded),
          thread_rank(thread_rank),
          thread_size(thread_size)
        {}
      };
    } // namespace detail
  } // namespace bcond
} // namespace libmpdataxx
