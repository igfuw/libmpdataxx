/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "../../blitz.hpp"

#include "../../bcond/cyclic_1d.hpp"
#include "../../bcond/cyclic_2d.hpp"
#include "../../bcond/cyclic_3d.hpp"
// TODO: split into 1D, 2D and 3D files?

#include "sharedmem.hpp"

namespace advoocat
{
  namespace concurr
  {
    template <typename real_t, int n_dims>
    class any
    {
      public:
      virtual void advance(int) { assert(false); }  
      virtual blitz::Array<real_t, n_dims> state(int e = 0) { assert(false); }
      virtual blitz::Array<real_t, n_dims> courant(int d = 0) { assert(false); }
    };

    namespace detail
    {
      template <
        class solver_t, 
        bcond::bcond_e bcx,
        bcond::bcond_e bcy,
        bcond::bcond_e bcz
      >
      class concurr_common : public any<typename solver_t::real_t, solver_t::n_dims>
      {
        // helper method to define subdomain ranges
	int min(int span, int rank, int size) 
	{ 
	  return rank * span / size; 
	}

	int max(int span, int rank, int size) 
	{ 
          return min(span, rank + 1, size) - 1; 
	}

	protected:

        // (cannot be nested due to tempaltes)
        typedef sharedmem<
          typename solver_t::real_t,
          solver_t::n_dims,
          solver_t::n_eqs,
          solver_t::n_tlev
        > mem_t;

	// member fields
	boost::ptr_vector<solver_t> algos; 
        std::unique_ptr<mem_t> mem;

	public:

        typedef typename solver_t::real_t real_t;

	// 1D ctor
	concurr_common(
	  const int s0, 
	  const typename solver_t::params_t &params,
          mem_t *mem_p,
	  const int n0
	)
	{
          mem.reset(mem_p);
	  solver_t::alloc(mem.get(), s0);

          std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr, shrd; // TODO: solver_t::bc_p

          if (bcx == bcond::cyclic)  // TODO: make a function taht does it
          {
            bxl.reset(new bcond::cyclic_left_1d<real_t>(rng_t(0, s0-1), solver_t::halo));
            bxr.reset(new bcond::cyclic_rght_1d<real_t>(rng_t(0, s0-1), solver_t::halo));
          }
          else assert(false);

	  for (int i0 = 0; i0 < n0; ++i0) 
          {
            shrd.reset(new bcond::shared<real_t>());
            const rng_t i(min(s0, i0, n0), max(s0, i0, n0)); 
	    algos.push_back(new solver_t(mem.get(), 
              i0 == 0      ? bxl : shrd,
              i0 == n0 - 1 ? bxr : shrd,
              i, params
            ));
          }
	}

	// 2D ctor
	concurr_common(
	  const int s0, const int s1, 
	  const typename solver_t::params_t &params,
          mem_t *mem_p,
	  const int n0, const int n1
	)
	{
          mem.reset(mem_p);
          solver_t::alloc(mem.get(), s0, s1);


// TODO: assert parallelisation in the right dimensions! (blitz::assertContiguous)
          for (int i0 = 0; i0 < n0; ++i0) 
          {
            for (int i1 = 0; i1 < n1; ++i1) 
            {
	      std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr, byl, byr, shrd;

	      if (bcx == bcond::cyclic)  // TODO: make a function taht does it
	      {
		bxl.reset(new bcond::cyclic_left_2d<0, real_t>(rng_t(0, s0-1), solver_t::halo));
		bxr.reset(new bcond::cyclic_rght_2d<0, real_t>(rng_t(0, s0-1), solver_t::halo));
	      } 
	      else assert(false);

	      if (bcy == bcond::cyclic)  // TODO: make a function taht does it
	      {
		byl.reset(new bcond::cyclic_left_2d<1, real_t>(rng_t(0, s1-1), solver_t::halo));
		byr.reset(new bcond::cyclic_rght_2d<1, real_t>(rng_t(0, s1-1), solver_t::halo));
	      }
	      else assert(false);

              shrd.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1

              const rng_t 
                i( min(s0, i0, n0), max(s0, i0, n0) ),
                j( min(s1, i1, n1), max(s1, i1, n1) );
              algos.push_back(new solver_t(mem.get(), 
                i0 == 0      ? bxl : shrd,
                i0 == n0 - 1 ? bxr : shrd,
                byl, 
                byr, 
                i, j, params
              ));
            }
          }
	}

	// 3D ctor
	concurr_common(
	  const int s0, const int s1, const int s2, 
	  const typename solver_t::params_t &params,
          mem_t *mem_p,
	  const int n0, const int n1, const int n2
	)
	{
          mem.reset(mem_p);
	  solver_t::alloc(mem.get(), s0, s1, s2);

	  std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr, byl, byr, bzl, bzr;

	  if (bcx == bcond::cyclic) // TODO: make a function that does it
	  {
	    bxl.reset(new bcond::cyclic_left_3d<0, real_t>(rng_t(0, s0-1), solver_t::halo));
	    bxr.reset(new bcond::cyclic_rght_3d<0, real_t>(rng_t(0, s0-1), solver_t::halo));
	  }
	  else assert(false);

	  if (bcy == bcond::cyclic) // TODO: make a function taht does it
	  {
	    byl.reset(new bcond::cyclic_left_3d<1, real_t>(rng_t(0, s1-1), solver_t::halo));
	    byr.reset(new bcond::cyclic_rght_3d<1, real_t>(rng_t(0, s1-1), solver_t::halo));
	  }
	  else assert(false);

	  if (bcz == bcond::cyclic)  // TODO: make a function taht does it
	  {
	    bzl.reset(new bcond::cyclic_left_3d<2, real_t>(rng_t(0, s2-1), solver_t::halo));
	    bzr.reset(new bcond::cyclic_rght_3d<2, real_t>(rng_t(0, s2-1), solver_t::halo));
	  }
	  else assert(false);

// TODO: parallel!
	  for (int i0 = 0; i0 < n0; ++i0) 
          {
	    for (int i1 = 0; i1 < n1; ++i1) 
            {
	      for (int i2 = 0; i2 < n2; ++i2) 
              {
                rng_t
                  i( min(s0, i0, n0), max(s0, i0, n0) ),
                  j( min(s1, i1, n1), max(s1, i1, n1) ),
                  k( min(s2, i2, n2), max(s2, i2, n2) );
		algos.push_back(new solver_t(mem.get(), bxl, bxr, byl, byr, bzl, bzr, i, j, k, params));
              }
            }
          }
        }

        virtual void solve(int nt) = 0;
    
        void advance(int nt) 
        {   
          solve(nt);
          mem->cycle();
        }  

	typename solver_t::arr_t state(int e = 0)
	{
	  return mem->state(e);
	}

	typename solver_t::arr_t courant(int d = 0)
	{
	  return mem->courant(d);
	}
      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace advoocat
