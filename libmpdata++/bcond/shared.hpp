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

      virtual void fill_halos_sclr(const blitz::Array<real_t, 1> &, const bool) { };
      virtual void fill_halos_sclr(const blitz::Array<real_t, 2> &, const rng_t &, const bool) { };
      virtual void fill_halos_sclr(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &, const bool) { };

      virtual void fill_halos_pres(const blitz::Array<real_t, 2> &, const rng_t &) { };
      virtual void fill_halos_pres(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) { };
      
      virtual void save_edge_vel(const blitz::Array<real_t, 2> &, const rng_t &) { };
      virtual void save_edge_vel(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) { };

      virtual void set_edge_pres(const blitz::Array<real_t, 2> &, const rng_t &, int) { };
      virtual void set_edge_pres(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &, int) { };

      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 1>> &) { };
      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 2>> &, const rng_t &) { };
      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 3>> &, const rng_t &, const rng_t &) { }; 

      virtual void fill_halos_vctr_nrml(const blitz::Array<real_t, 2> &, const rng_t &) { };
      virtual void fill_halos_vctr_nrml(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) { };
    };
  } // namespace bcond
} // namespace libmpdataxx
