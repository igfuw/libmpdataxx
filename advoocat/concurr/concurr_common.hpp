/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "../blitz.hpp"

#include "../bcond/cyclic_1d.hpp"
#include "../bcond/cyclic_2d.hpp"
#include "../bcond/cyclic_3d.hpp"
// TODO: split into 1D, 2D and 3D files?

namespace advoocat
{
  namespace concurr
  {
    namespace detail
    {
      template <
        class solver_t, 
        bcond::bcond_e bcx,
        bcond::bcond_e bcy,
        bcond::bcond_e bcz
      >
      class concurr_common
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

	// member fields
	typename solver_t::mem_t mem;
       
	protected:

        // nested class
        struct mem_t {};

	// member fields
	boost::ptr_vector<solver_t> algos; 

	public:

	typedef typename solver_t::mem_t::real_t real_t; 

	// 1D ctor
	concurr_common(
	  const int s0, 
	  const typename solver_t::params_t &params,
	  const int n0
	)
	  : mem(s0)
	{
	  solver_t::alloc(mem, s0);

          std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr;

          if (bcx == bcond::cyclic)  // TODO: make a function taht does it
          {
            bxl.reset(new bcond::cyclic_left_1d<real_t>(rng_t(0, s0-1), solver_t::halo));
            bxr.reset(new bcond::cyclic_rght_1d<real_t>(rng_t(0, s0-1), solver_t::halo));
          }
          else assert(false);

	  for (int i0 = 0; i0 < n0; ++i0) 
          {
            const rng_t i(min(s0, i0, n0), max(s0, i0, n0)); 
	    algos.push_back(new solver_t(mem, bxl, bxr, i, params));
          }
	}

	// 2D ctor
	concurr_common(
	  const int s0, const int s1, 
	  const typename solver_t::params_t &params,
	  const int n0, const int n1
	)
	  : mem(s0, s1)
	{
          solver_t::alloc(mem, s0, s1);

	  std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr, byl, byr;

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

          for (int i0 = 0; i0 < n0; ++i0) 
          {
            for (int i1 = 0; i1 < n1; ++i1) 
            {
              const rng_t 
                i( min(s0, i0, n0), max(s0, i0, n0) ),
                j( min(s1, i1, n1), max(s1, i1, n1) );
              algos.push_back(new solver_t(mem, bxl, bxr, byl, byr, i, j, params));
            }
          }
	}

	// 3D ctor
	concurr_common(
	  const int s0, const int s1, const int s2, 
	  const typename solver_t::params_t &params,
	  const int n0, const int n1, const int n2
	)
	  : mem(s0, s1, s2)
	{
	  solver_t::alloc(mem, s0, s1, s2);

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
		algos.push_back(new solver_t( mem, bxl, bxr, byl, byr, bzl, bzr, i, j, k, params));
              }
            }
          }
        }

	typename solver_t::mem_t::arr_t state(int e = 0)
	{
	  return mem.state(e);
	}

	typename solver_t::mem_t::arr_t courant(int d = 0)
	{
	  return mem.courant(d);
	}
      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace advoocat
