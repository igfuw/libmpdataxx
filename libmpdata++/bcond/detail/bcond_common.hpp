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

    enum bcond_e { null, cyclic, polar, open, rigid, remote }; 
    enum drctn_e { left, rght };

    template<
      typename real_t, 
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
      template <typename real_t>
      class bcond_common
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
	
        virtual void fill_halos_pres(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &)
	{
	  assert(false && "bcond::fill_halos_pres() called!");
	};
	
	// setting to a value
	virtual void set_edge_pres(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) 
	{
	  assert(false && "bcond::set_edge_a() called!");
	};
	
	// setting to a given array
	virtual void set_edge_pres(const blitz::Array<real_t, 3> &,
                                   const blitz::Array<real_t, 3> &,
                                   const rng_t &,
                                   const rng_t &
                                  ) 
	{
	  assert(false && "bcond::set_edge_b() called!");
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
	bcond_common(const rng_t &i, const int halo, const int grid_size_0) :
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
    }; // namespace detail
  }; // namespace bcond
}; // namespace libmpdataxx
