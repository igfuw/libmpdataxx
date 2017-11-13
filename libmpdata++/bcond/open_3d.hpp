// 3D open boundary conditions for libmpdata++
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
        n_dims == 3
      >::type 
    > : public detail::bcond_common<real_t, halo>
    { 
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor
      
      // holds saved initial value of edge velocity
      arr_t edge_velocity;

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
	using namespace idxperm;
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
        {
          if (deriv)
	    a(pi<d>(i, j, k)) = 0;
          else
	    a(pi<d>(i, j, k)) = a(pi<d>(this->left_edge_sclr, j, k));
        }
      }
      
      void fill_halos_pres(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        a(pi<d>(this->left_halo_sclr.last(), j, k)) = 2 * a(pi<d>(this->left_edge_sclr,     j, k))
                                                        - a(pi<d>(this->left_edge_sclr + 1, j, k));
        if (halo > 1)
        {
          a(pi<d>(this->left_halo_sclr.last() - 1, j, k)) =   3 * a(pi<d>(this->left_edge_sclr,     j, k))
                                                            - 2 * a(pi<d>(this->left_edge_sclr + 1, j, k));
        }
      }
      
      void save_edge_vel(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        auto s = a.shape();
        s[d] = 1;
        edge_velocity.resize(s);
        edge_velocity(pi<d>(0, j, k)) = a(pi<d>(this->left_edge_sclr, j, k));
      }
      
      void set_edge_pres(arr_t &a, const rng_t &j, const rng_t &k, int sign)
      {
        using namespace idxperm;
        a(pi<d>(this->left_edge_sclr, j, k)) = sign * edge_velocity(pi<d>(0, j, k));
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
	using namespace idxperm;
	const int i = this->left_edge_sclr;

        // TODO: exactly the same code below!
        switch (d) // note: order and lack of breaks intentional!
        {
          case 0:
          av[d+2](pi<d>(i, j, (k-h).first())) = 0;
          av[d+2](pi<d>(i, j, (k+h).last() )) = 0;

          case 1:
          av[d+1](pi<d>(i, (j-h).first(), k)) = 0;
          av[d+1](pi<d>(i, (j+h).last(),  k)) = 0;

          case 2: 
          break;

          default: assert(false);
        }

	assert(std::isfinite(sum(av[d  ](pi<d>(i+h, j, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j-h, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j+h, k)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k-h)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k+h)))));

        // zero-divergence condition
        for (int ii = this->left_halo_vctr.first(); ii <= this->left_halo_vctr.last() - (ad ? 1 : 0); ++ii)
        {
          av[d](pi<d>(ii, j, k)) = 
            av[d](pi<d>(i+h, j, k)) 
            -(
              av[d+1](pi<d>(i, j-h, k)) - 
              av[d+1](pi<d>(i, j+h, k))   
            ) 
            -(
              av[d+2](pi<d>(i, j, k-h)) -
              av[d+2](pi<d>(i, j, k+h)) 
            );
        }
      }

      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
          a(pi<d>(i, j, k)) = 0; 
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == open &&
        dir == rght &&
        n_dims == 3
      >::type 
    > : public detail::bcond_common<real_t, halo>
    { 
      using parent_t = detail::bcond_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor
      
      // holds saved initial value of edge velocity
      arr_t edge_velocity;
      
      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
	using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
        {
          if (deriv)
	    a(pi<d>(i, j, k)) = 0;
          else
	    a(pi<d>(i, j, k)) = a(pi<d>(this->rght_edge_sclr, j, k));
        }
      }
      
      void fill_halos_pres(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // equivalent to one-sided derivatives at the boundary
        a(pi<d>(this->rght_halo_sclr.first(), j, k)) = 2 * a(pi<d>(this->rght_edge_sclr,     j, k))
                                                         - a(pi<d>(this->rght_edge_sclr - 1, j, k));

        if (halo > 1)
        {
          a(pi<d>(this->rght_halo_sclr.first() + 1, j, k)) =   3 * a(pi<d>(this->rght_edge_sclr,     j, k))
                                                             - 2 * a(pi<d>(this->rght_edge_sclr - 1, j, k));
        }
      }
      
      void save_edge_vel(const arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        auto s = a.shape();
        s[d] = 1;
        edge_velocity.resize(s);
        edge_velocity(pi<d>(0, j, k)) = a(pi<d>(this->rght_edge_sclr, j, k));
      }
      
      void set_edge_pres(arr_t &a, const rng_t &j, const rng_t &k, int sign)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_edge_sclr, j, k)) = sign * edge_velocity(pi<d>(0, j, k));
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k, const bool ad = false)
      {
	using namespace idxperm;
        const int i = this->rght_edge_sclr;

        switch (d) // note: order and lack of breaks intentional!
        {
          case 0:
	  av[d+2](pi<d>(i, j, (k-h).first())) = 0;
	  av[d+2](pi<d>(i, j, (k+h).last() )) = 0;

          case 1:
	  av[d+1](pi<d>(i, (j-h).first(), k)) = 0;
	  av[d+1](pi<d>(i, (j+h).last(),  k)) = 0;

          case 2: 
          break;

          default: assert(false);
        }

	assert(std::isfinite(sum(av[d  ](pi<d>(i-h, j, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j-h, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j+h, k)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k-h)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k+h)))));

        for (int ii = this->rght_halo_vctr.first() + (ad ? 1 : 0); ii <= this->rght_halo_vctr.last(); ++ii)
        {
          av[d](pi<d>(ii, j, k)) = 
            av[d](pi<d>(i-h, j, k)) 
            +(
              av[d+1](pi<d>(i, j-h, k)) - 
              av[d+1](pi<d>(i, j+h, k)) 
            )
            +(
              av[d+2](pi<d>(i, j, k-h)) - 
              av[d+2](pi<d>(i, j, k+h)) 
            );
        }
      }
      
      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
          a(pi<d>(i, j, k)) = 0; 
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
