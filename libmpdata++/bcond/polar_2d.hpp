/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv2+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/bcond/polar_common.hpp>
#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template<int d, typename real_t>
    class polar_left_2d : public polar_common<real_t>
    {
      using parent_t = polar_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      void bcinit(const arr_t &a, const rng_t &j) {}

      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = 0; i < this->halo; ++i)
        {
          a(pi<d>(this->left_halo_sclr.last() - i,
                  this->west_polar_halo(j)        ))
          =
          a(pi<d>(this->left_edge_sclr + i,
                  this->east_sclr_polar_edge(j) ));

          a(pi<d>(this->left_halo_sclr.last() - i,
                  this->east_polar_halo(j)        ))
          =
          a(pi<d>(this->left_edge_sclr + i,
                  this->west_sclr_polar_edge(j) ));
        }
      }

      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
	a(pi<d>(this->left_halo_vctr, j)) = 0;
      }

      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = 0; i < this->halo; ++i)
        {
          a(pi<d>(this->left_halo_sclr.last() - i,
                  this->east_polar_halo(j)        ))
          =
          a(pi<d>(this->left_intr_vctr.first() + i,
                     this->west_vctr_polar_edge(j)    ));
          
          a(pi<d>(this->left_halo_sclr.last() - i,
                  this->west_polar_halo(j)        ))
          =
          a(pi<d>(this->left_intr_vctr.first() + i,
                  this->east_vctr_polar_edge(j)    ));
        }
      }
    };

    template<int d, typename real_t>
    class polar_rght_2d : public polar_common<real_t>
    {
      using parent_t = polar_common<real_t>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // method invoked by the solver
      
      void bcinit(const arr_t &a, const rng_t &j) {}

      void fill_halos_sclr(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;

        for (int i = 0; i < this->halo; ++i)
        { 
          a(pi<d>(this->rght_halo_sclr.first() + i,
                  this->west_polar_halo(j)         ))
          =
          a(pi<d>(this->rght_edge_sclr - i,
                  this->east_sclr_polar_edge(j) ));

          a(pi<d>(this->rght_halo_sclr.first() + i,
                  this->east_polar_halo(j)         ))
          =
          a(pi<d>(this->rght_edge_sclr - i,
                  this->west_sclr_polar_edge(j) ));
        }
      }

      void fill_halos_vctr_alng(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
	a(pi<d>(this->rght_halo_vctr, j)) = 0;
      }
      
      void fill_halos_vctr_nrml(const arr_t &a, const rng_t &j)
      {
	using namespace idxperm;
        for (int i = 0; i < this->halo; ++i)
        {
          a(pi<d>(this->rght_halo_sclr.first() + i,
                  this->east_polar_halo(j)         ))
          =
          a(pi<d>(this->rght_intr_vctr.last() - i,
                  this->west_vctr_polar_edge(j)   ));
          
          a(pi<d>(this->rght_halo_sclr.first() + i,
                  this->west_polar_halo(j)         ))
          =
          a(pi<d>(this->rght_intr_vctr.last() - i,
                  this->east_vctr_polar_edge(j)   ));
        }
      }
    };
  }; // namespace bcond
}; // namespace libmpdataxx
