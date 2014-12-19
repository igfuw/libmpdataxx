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
    template<int d, typename real_t>
    class polar_left_2d : public detail::polar_common<real_t>
    {
      using parent_t = detail::polar_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void fill_halos_sclr(const arr_t &a, const rng_t &j, const bool deriv = false)
      {
	using namespace idxperm;
        for (int i = 0; i < this->halo; ++i)
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

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
	av[d](pi<d>(this->left_halo_vctr, j)) = 0;
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = 0; i < this->halo; ++i)
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

    template<int d, typename real_t>
    class polar_rght_2d : public detail::polar_common<real_t>
    {
      using parent_t = detail::polar_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void fill_halos_sclr(const arr_t &a, const rng_t &j, const bool deriv = false)
      {
	using namespace idxperm;

        for (int i = 0; i < this->halo; ++i)
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

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j)
      {
	using namespace idxperm;
	av[d](pi<d>(this->rght_halo_vctr, j)) = 0;
      }
      
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = 0; i < this->halo; ++i)
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
  }; // namespace bcond
}; // namespace libmpdataxx
