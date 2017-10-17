// do-nothing shared-memory boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo>
    class shared : public detail::bcond_common<real_t, halo>
    {
      public:
      
      using parent_t = detail::bcond_common<real_t, halo>;
      
      using arr_1d_t = typename parent_t::arr_1d_t;
      using arr_2d_t = typename parent_t::arr_2d_t;
      using arr_3d_t = typename parent_t::arr_3d_t;

      virtual void fill_halos_sclr(arr_1d_t &, const bool) { };
      virtual void fill_halos_sclr(arr_2d_t &, const rng_t &, const bool) { };
      virtual void fill_halos_sclr(arr_3d_t &, const rng_t &, const rng_t &, const bool) { };

      virtual void fill_halos_pres(arr_2d_t &, const rng_t &) { };
      virtual void fill_halos_pres(arr_3d_t &, const rng_t &, const rng_t &) { };
      
      virtual void save_edge_vel(const arr_2d_t &, const rng_t &) { };
      virtual void save_edge_vel(const arr_3d_t &, const rng_t &, const rng_t &) { };

      virtual void set_edge_pres(arr_2d_t &, const rng_t &, int) { };
      virtual void set_edge_pres(arr_3d_t &, const rng_t &, const rng_t &, int) { };

      virtual void fill_halos_vctr_alng(arrvec_t<arr_1d_t> &, const bool) { };
      virtual void fill_halos_vctr_alng(arrvec_t<arr_2d_t> &, const rng_t &, const bool) { };
      virtual void fill_halos_vctr_alng(arrvec_t<arr_3d_t> &, const rng_t &, const rng_t &, const bool) { }; 
      
      virtual void fill_halos_vctr_nrml(arr_2d_t &, const rng_t &) { };
      virtual void fill_halos_vctr_nrml(arr_3d_t &, const rng_t &, const rng_t &) { };
      
      virtual void fill_halos_sgs_div(arr_2d_t &, const rng_t &) { };
      virtual void fill_halos_sgs_div(arr_3d_t &, const rng_t &, const rng_t &) { };
      
      virtual void fill_halos_sgs_vctr(arrvec_t<arr_2d_t> &, const arr_2d_t &, const rng_t &, const int offset = 0) { };
      virtual void fill_halos_sgs_vctr(arrvec_t<arr_3d_t> &, const arr_3d_t &, const rng_t &, const rng_t &, const int offset = 0) { };
      
      virtual void fill_halos_sgs_tnsr(arrvec_t<arr_2d_t> &, const arr_2d_t &, const arr_2d_t &, const rng_t &, const real_t) { };
      virtual void fill_halos_sgs_tnsr(arrvec_t<arr_3d_t> &, const arr_3d_t &, const arr_3d_t &,
                                                             const rng_t &, const rng_t &, const real_t) { };

      using parent_t::parent_t;
      // ctor
      // parent constructor takes minimal parameters that construct valid member ranges, note that they
      // aren't actually used !
      shared() : parent_t(rng_t(0, 2), 2) {} // TODO: bcond_any?
    };
  } // namespace bcond
} // namespace libmpdataxx
