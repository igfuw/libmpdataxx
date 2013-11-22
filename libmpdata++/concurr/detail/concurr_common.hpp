/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include <libmpdata++/blitz.hpp>

#include <libmpdata++/bcond/cyclic_1d.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/bcond/cyclic_3d.hpp>
// TODO: split into 1D, 2D and 3D files?

#include <libmpdata++/concurr/detail/sharedmem.hpp>

namespace libmpdataxx
{
  /// @brief concurr namespace
  namespace concurr
  {
    template <typename real_t_, int n_dims>
    class any
    {
      public:
      typedef real_t_ real_t;
      virtual void advance(int) { assert(false); throw; }  
      virtual blitz::Array<real_t_, n_dims> advectee(int e = 0) { assert(false); throw; }
      virtual blitz::Array<real_t_, n_dims> advector(int d = 0) { assert(false); throw; }
      virtual blitz::Array<real_t_, n_dims> g_factor() { assert(false); throw; }
      virtual bool *panic_ptr() { assert(false && "unimplemented!"); throw; }
    };

    namespace detail
    {
      template <
        class solver_t_, 
        bcond::bcond_e bcx,
        bcond::bcond_e bcy,
        bcond::bcond_e bcz
      >
      class concurr_common : public any<typename solver_t_::real_t, solver_t_::n_dims>
      {
        public:

        typedef solver_t_ solver_t;

        private:

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
          solver_t::n_tlev
        > mem_t;

	// member fields
	boost::ptr_vector<solver_t> algos; 
        std::unique_ptr<mem_t> mem;

	public:

        typedef typename solver_t::real_t real_t;

	// ctor
	concurr_common(
	  const typename solver_t::params_t &p,
          mem_t *mem_p,
	  const int &size
	) {
          mem.reset(mem_p);
	  solver_t::alloc(mem.get(), p);
          init(p, p.span, size); 
        }

        private:
 
        void init(
          const typename solver_t::params_t &p,
          const std::array<int, 1> &span, const int &n0
        )
        {
          std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr, shrdl, shrdr; // TODO: solver_t::bc_p

          if (bcx == bcond::cyclic)  // TODO: make a function that does it
          {
            bxl.reset(new bcond::cyclic_left_1d<real_t>(rng_t(0, span[0]-1), solver_t::halo));
            bxr.reset(new bcond::cyclic_rght_1d<real_t>(rng_t(0, span[0]-1), solver_t::halo));
          }
          else assert(false);

	  for (int i0 = 0; i0 < n0; ++i0) 
          {
            shrdl.reset(new bcond::shared<real_t>());
            shrdr.reset(new bcond::shared<real_t>());
            const rng_t i(min(span[0], i0, n0), max(span[0], i0, n0)); 
	    algos.push_back(
              new solver_t(
                typename solver_t::ctor_args_t({
		  mem.get(), 
		  i0 == 0      ? bxl : shrdl,
		  i0 == n0 - 1 ? bxr : shrdr,
		  i
                }), 
                p
              )
            );
          }
	}

        void init(
          const typename solver_t::params_t &p,
	  const std::array<int, 2> &span, 
          const int &n0, const int &n1 = 1
        ) {
// TODO: assert parallelisation in the right dimensions! (blitz::assertContiguous)
          for (int i0 = 0; i0 < n0; ++i0) 
          {
            for (int i1 = 0; i1 < n1; ++i1) 
            {
	      std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr, byl, byr, shrdl, shrdr;

	      if (bcx == bcond::cyclic)  // TODO: make a function taht does it
	      {
		bxl.reset(new bcond::cyclic_left_2d<0, real_t>(rng_t(0, span[0]-1), solver_t::halo));
		bxr.reset(new bcond::cyclic_rght_2d<0, real_t>(rng_t(0, span[0]-1), solver_t::halo));
	      } 
	      else assert(false);

	      if (bcy == bcond::cyclic)  // TODO: make a function taht does it
	      {
		byl.reset(new bcond::cyclic_left_2d<1, real_t>(rng_t(0, span[1]-1), solver_t::halo));
		byr.reset(new bcond::cyclic_rght_2d<1, real_t>(rng_t(0, span[1]-1), solver_t::halo));
	      }
	      else assert(false);

              shrdl.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1
              shrdr.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1

              const rng_t 
                i( min(span[0], i0, n0), max(span[0], i0, n0) ),
                j( min(span[1], i1, n1), max(span[1], i1, n1) );
              algos.push_back(
                new solver_t(
                  typename solver_t::ctor_args_t({
		    mem.get(), 
		    i0 == 0      ? bxl : shrdl,
		    i0 == n0 - 1 ? bxr : shrdr,
		    byl, 
		    byr, 
		    i, j
                  }), 
                  p
                )
              );
            }
          }
	}

        void init(
          const typename solver_t::params_t &p,
	  const std::array<int, 3> &span, 
          const int &n0, const int &n1 = 1, const int &n2 = 1
        ) {
	  std::unique_ptr<bcond::bcond_t<real_t>> bxl, bxr, byl, byr, bzl, bzr, shrdl, shrdr;

// TODO: renew pointers only if invalid ?
	  for (int i0 = 0; i0 < n0; ++i0) 
          {
	    for (int i1 = 0; i1 < n1; ++i1) 
            {
	      for (int i2 = 0; i2 < n2; ++i2) 
              {
                if (bcx == bcond::cyclic) // TODO: make a function that does it
                {
                  bxl.reset(new bcond::cyclic_left_3d<0, real_t>(rng_t(0, span[0]-1), solver_t::halo));
                  bxr.reset(new bcond::cyclic_rght_3d<0, real_t>(rng_t(0, span[0]-1), solver_t::halo));
                }
                else assert(false);

                if (bcy == bcond::cyclic) // TODO: make a function taht does it
                {
                  byl.reset(new bcond::cyclic_left_3d<1, real_t>(rng_t(0, span[1]-1), solver_t::halo));
                  byr.reset(new bcond::cyclic_rght_3d<1, real_t>(rng_t(0, span[1]-1), solver_t::halo));
                }
                else assert(false);

                if (bcz == bcond::cyclic)  // TODO: make a function taht does it
                {
                  bzl.reset(new bcond::cyclic_left_3d<2, real_t>(rng_t(0, span[2]-1), solver_t::halo));
                  bzr.reset(new bcond::cyclic_rght_3d<2, real_t>(rng_t(0, span[2]-1), solver_t::halo));
                }
                else assert(false);

                shrdl.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1
                shrdr.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1

                rng_t
                  i( min(span[0], i0, n0), max(span[0], i0, n0) ),
                  j( min(span[1], i1, n1), max(span[1], i1, n1) ),
                  k( min(span[2], i2, n2), max(span[2], i2, n2) );

		algos.push_back(
                  new solver_t(
                    typename solver_t::ctor_args_t({
                      mem.get(), 
		      i0 == 0      ? bxl : shrdl,
		      i0 == n0 - 1 ? bxr : shrdr,
                      byl, byr, 
                      bzl, bzr, 
                      i, j, k
                    }), 
                    p
                  )
                );
              }
            }
          }
        }

        public:

        virtual void solve(int nt) = 0;
    
        void advance(int nt) 
        {   
          solve(nt);
        }  

	typename solver_t::arr_t advectee(int e = 0)
	{
	  return mem->advectee(e);
	}

	typename solver_t::arr_t advector(int d = 0)
	{
	  return mem->advector(d);
	}

	typename solver_t::arr_t g_factor()
	{
	  return mem->g_factor();
	}

        bool *panic_ptr()
        {
          return &this->mem->panic;
        }
      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace libmpdataxx
