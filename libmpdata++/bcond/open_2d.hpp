// 2D open boundary conditions for libmpdata++
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
        knd == open &&
        dir == left &&
        n_dims == 2
      >::type
    > : public detail::bcond_common<real_t, halo>
    {
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const bool deriv = false)
      {
	using namespace idxperm;
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
        {
          if (deriv)
	    a(pi<d>(i, j)) = 0;
          else 
	    a(pi<d>(i, j)) = a(pi<d>(this->left_edge_sclr, j)); // zero-gradient condition for scalar
        }
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
	const int i = this->left_edge_sclr;
   
        // if executed first (d=0) this could contain NaNs
        if (d == 0) 
        {
          av[d+1](pi<d>(i, (j-h).first())) = 0;
          av[d+1](pi<d>(i, (j+h).last())) = 0;
        }
       
	// zero-divergence condition
        for (int ii = this->left_halo_vctr.first(); ii <= this->left_halo_vctr.last(); ++ii)
        {
	  av[d](pi<d>(ii, j)) = 
	    av[d](pi<d>(i+h, j)) 
            -(
	      av[d+1](pi<d>(i, j-h)) -
	      av[d+1](pi<d>(i, j+h))
	    );
        }
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
          a(pi<d>(i, j)) = 0; 
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == open &&
        dir == rght &&
        n_dims == 2
      >::type
    > : public detail::bcond_common<real_t, halo>
    {
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor
      
      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const bool deriv = false)
      {
	using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
        {
	  if (deriv)
            a(pi<d>(i, j)) = 0; // zero gradient for scalar gradient
          else
            a(pi<d>(i, j)) = a(pi<d>(this->rght_edge_sclr, j)); // zero gradient for scalar
        }
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
	const int i = this->rght_edge_sclr;

        // if executed first (d=0) this could contain NaNs
        if (d == 0) 
        {
          av[d+1](pi<d>(i, (j-h).first())) = 0;
          av[d+1](pi<d>(i, (j+h).last())) = 0;
        }
       
	// zero-divergence condition
        for (int ii = this->rght_halo_vctr.first(); ii <= this->rght_halo_vctr.last(); ++ii)
        {
	  av[d](pi<d>(ii, j)) = 
	    av[d](pi<d>(i-h, j)) + (
	      av[d+1](pi<d>(i, j-h)) -
	      av[d+1](pi<d>(i, j+h))
	    );
        }
      }
      
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
          a(pi<d>(i, j)) = 0; 
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
