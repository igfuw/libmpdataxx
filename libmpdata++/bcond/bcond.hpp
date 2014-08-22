/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  */

#pragma once

#include <libmpdata++/formulae/arakawa_c.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    using namespace arakawa_c;

    enum bcond_e { null, cyclic, polar, open, rigid }; 

    template <typename real_t>
    class bcond_t
    {
      public:


      // 1D
      virtual void fill_halos_sclr(const blitz::Array<real_t, 1> &, const bool deriv = false) 
      { 
        assert(false && "bcond::fill_halos_sclr() called!"); 
      };
      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 1>> &) 
      { 
        assert(false && "bcond::fill_halos_vctr() called!"); 
      };

      // 2D
      virtual void fill_halos_sclr(const blitz::Array<real_t, 2> &, const rng_t &, const bool deriv = false) 
      {
        assert(false && "bcond::fill_halos_sclr() called!");
      };
      
      virtual void fill_halos_pres(const blitz::Array<real_t, 2> &, const rng_t &)
      {
        assert(false && "bcond::fill_halos_pres() called!");
      };
      
      // setting to a value
      virtual void set_edge_pres(const blitz::Array<real_t, 2> &, const rng_t &) 
      {
        assert(false && "bcond::set_edge_a() called!");
      };
      
      // setting to a given array
      virtual void set_edge_pres(const blitz::Array<real_t, 2> &,const blitz::Array<real_t, 2> &, const rng_t &) 
      {
        assert(false && "bcond::set_edge_b() called!");
      };

      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 2>> &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_vctr_alng() called!");
      };
      virtual void fill_halos_vctr_nrml(const blitz::Array<real_t, 2> &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_vctr_nrml() called!");
      };

      // 3D
      virtual void fill_halos_sclr(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &, const bool deriv = false) 
      {
        assert(false && "bcond::fill_halos_sclr() called!");
      };
      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 3>> &, const rng_t &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_vctr() called!");
      };
      virtual void fill_halos_vctr_nrml(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_vctr_nrml() called!");
      };
      

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
      const int halo;

      public:

      // ctor
      bcond_t(const rng_t &i, const int halo) :
        halo(halo),
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
        )
      {} 
    };

    template <typename real_t>
    class shared : public bcond_t<real_t> // TODO: move to a bcond_shared file and document!
    {
      public:

      virtual void fill_halos_sclr(const blitz::Array<real_t, 1> &, const bool) { };
      virtual void fill_halos_sclr(const blitz::Array<real_t, 2> &, const rng_t &, const bool) { };
      virtual void fill_halos_sclr(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &, const bool) { };

      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 1>> &) { };
      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 2>> &, const rng_t &) { };
      virtual void fill_halos_vctr_alng(const arrvec_t<blitz::Array<real_t, 3>> &, const rng_t &, const rng_t &) { }; 

      virtual void fill_halos_vctr_nrml(const blitz::Array<real_t, 2> &, const rng_t &) { };
      virtual void fill_halos_vctr_nrml(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) { };
      
      using parent_t = bcond_t<real_t>;
      using parent_t::parent_t;

      // ctor
      // parent constructor takes minimal parameters that construct valid member ranges, note that they
      // aren't actually used !
      shared() : parent_t(rng_t(0, 2), 1) {}
    };

  }; // namespace bcond
}; // namespace libmpdataxx
