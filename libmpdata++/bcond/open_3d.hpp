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
    class open_left_3d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
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

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
	const int i = this->left_edge_sclr;

        switch (d) // note: order and lack of breaks intentional!
        {
          case 0:
          av[d+2](pi<d>(i, j, (k-h).first())) = 0;
          av[d+2](pi<d>(i, j, (k+h).last() )) = 0;

          case 1:
          av[d+1](pi<d>(i, (j-h).first(), k)) = 0;
          av[d+1](pi<d>(i, (j+h).last(),  k)) = 0;

          default: assert(false);
        }

	assert(std::isfinite(sum(av[d  ](pi<d>(i+h, j, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j-h, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j+h, k)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k-h)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k+h)))));

        // zero-divergence condition
        for (int ii = this->left_halo_vctr.first(); ii <= this->left_halo_vctr.last(); ++ii)
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

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
          a(pi<d>(i, j, k)) = 0; 
      }
    };

    template<int d, typename real_t>
    class open_rght_3d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor
      
      public:

      void fill_halos_sclr(const arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
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

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av, const rng_t &j, const rng_t &k)
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

          default: assert(false);
        }

	assert(std::isfinite(sum(av[d  ](pi<d>(i-h, j, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j-h, k)))));
	assert(std::isfinite(sum(av[d+1](pi<d>(i, j+h, k)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k-h)))));
	assert(std::isfinite(sum(av[d+2](pi<d>(i, j, k+h)))));

        for (int ii = this->rght_halo_vctr.first(); ii <= this->rght_halo_vctr.last(); ++ii)
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
      
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j, const rng_t &k)
      {
	using namespace idxperm;
        // note intentional sclr
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
          a(pi<d>(i, j, k)) = 0; 
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
