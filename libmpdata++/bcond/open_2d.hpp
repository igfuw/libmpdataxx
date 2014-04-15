/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/bcond/open_common.hpp>
#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template<int d, typename real_t>
    class open_left_2d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
	  a(pi<d>(i, j)) = a(pi<d>(this->left_edge_sclr, j));
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

    template<int d, typename real_t>
    class open_rght_2d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor
      
      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
	  a(pi<d>(i, j)) = a(pi<d>(this->rght_edge_sclr, j));
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
