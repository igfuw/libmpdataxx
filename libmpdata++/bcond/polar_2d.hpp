// 2D polar boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/polar_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == polar &&
        dir == left &&
        n_dims == 2
      >::type 
    > : public detail::polar_common<real_t, halo>
    { 
      using parent_t = detail::polar_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
	using namespace idxperm;
        for (int i = 0; i < halo; ++i)
        { 
          for (int jj = j.first(); jj <= j.last(); jj++)
          {
            a(pi<d>(this->left_halo_sclr.last() - i,
                    jj)) 
            =
            a(pi<d>(this->left_edge_sclr + i,
                    this->polar_neighbours(jj)));
              
          }
        }
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
	av[d](pi<d>(this->left_halo_vctr.last(), j)) = 0;
        if (halo > 1)
        {
          for (int jj = j.first(); jj <= j.last(); jj++)
          {
	    av[d](pi<d>(this->left_halo_vctr.first(), jj))
            = 
            av[d](pi<d>(this->left_edge_sclr + h, this->polar_neighbours(jj)));
          }
        }
      }

      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = 0; i < halo; ++i)
        { 
          for (int jj = j.first(); jj <= j.last(); jj++)
          {
            a(pi<d>(this->left_halo_sclr.first() + i,
                    jj + h)) 
            =
            a(pi<d>(this->left_intr_vctr.last() - i,
                    this->polar_neighbours(jj) + h));
              
          }
        }
      }
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == polar &&
        dir == rght &&
        n_dims == 2
      >::type
    > : public detail::polar_common<real_t, halo>
    { 
      using parent_t = detail::polar_common<real_t, halo>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
	using namespace idxperm;

        for (int i = 0; i < halo; ++i)
        { 
          for (int jj = j.first(); jj <= j.last(); jj++)
          {
            a(pi<d>(this->rght_halo_sclr.first() + i,
                    jj)) 
            =
            a(pi<d>(this->rght_edge_sclr - i,
                    this->polar_neighbours(jj)));
              
          }
        }
      }

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
	av[d](pi<d>(this->rght_halo_vctr.first(), j)) = 0;
        if (halo > 1)
        {
          for (int jj = j.first(); jj <= j.last(); jj++)
          {
	    av[d](pi<d>(this->rght_halo_vctr.last(), jj))
            =
	    av[d](pi<d>(this->rght_edge_sclr - h, this->polar_neighbours(jj)));
          }
        }
      }
      
      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = 0; i < halo; ++i)
        { 
          for (int jj = j.first(); jj <= j.last(); jj++)
          {
            a(pi<d>(this->rght_halo_sclr.first() + i,
                    jj + h)) 
            =
            a(pi<d>(this->rght_intr_vctr.last() - i,
                    this->polar_neighbours(jj) + h));
              
          }
        }
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
