/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include <libmpdata++/blitz.hpp>

#include <libmpdata++/bcond/shared.hpp>
#include <libmpdata++/bcond/cyclic_1d.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/bcond/cyclic_3d.hpp>
#include <libmpdata++/bcond/open_1d.hpp>
#include <libmpdata++/bcond/open_2d.hpp>
#include <libmpdata++/bcond/open_3d.hpp>
#include <libmpdata++/bcond/polar_2d.hpp>
#include <libmpdata++/bcond/rigid_2d.hpp>

#include <libmpdata++/concurr/detail/sharedmem.hpp>
#include <libmpdata++/concurr/detail/timer.hpp>

// TODO: simplify 1d/2d/3d logic below or split into separate files?

namespace libmpdataxx
{
  /// @brief concurr namespace
  namespace concurr
  {
    template <typename real_t, int n_dims>
    struct any
    {
      virtual 
//<listing-1>
      void advance(int) 
//</listing-1>
      { assert(false); throw; }  

      virtual 
//<listing-2>
      blitz::Array<real_t, n_dims> advectee(int eqn = 0)
//</listing-2>
      { assert(false); throw; }

      virtual 
//<listing-3>
      blitz::Array<real_t, n_dims> advector(int dim = 0) 
//</listing-3>
      { assert(false); throw; }

      virtual 
//<listing-4>
      blitz::Array<real_t, n_dims> g_factor() 
//</listing-4>
      { assert(false); throw; }

      virtual 
//<listing-5>
      bool *panic_ptr() 
//</listing-5>
      { assert(false && "unimplemented!"); throw; }

      // dtor
      virtual ~any() {}
    };

    namespace detail
    {
      template<
        class solver_t_, 
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
      >
      class concurr_common : public any<typename solver_t_::real_t, solver_t_::n_dims>
      {

        public:

        typedef solver_t_ solver_t;
        
        static_assert(
          (solver_t::n_dims == 3) ||
          (solver_t::n_dims == 2 
            && bczl == bcond::null 
            && bczr == bcond::null
          ) ||
          (solver_t::n_dims == 1 
            && bczl == bcond::null 
            && bczr == bcond::null
            && bcyl == bcond::null 
            && bcyr == bcond::null
          )
          ,
          "more boundary conditions than dimensions"
        );

        private:

        // helper methods to define subdomain ranges
	int min(const int &grid_size, const int &rank, const int &size) 
	{ 
	  return rank * grid_size / size; 
	}

	int max(const int &grid_size, const int &rank, const int &size) 
	{ 
          return min(grid_size, rank + 1, size) - 1; 
	}

        template <int d>
        rng_t slab(
          const std::array<int, solver_t::n_dims> &grid_size, 
          const int &rank = 0, 
          const int &size = 1
        ) {
          return rng_t(
            min(grid_size[d], rank, size),
            max(grid_size[d], rank, size)
          );
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
        timer tmr;

	public:

        typedef typename solver_t::real_t real_t;

        // dtor
	virtual ~concurr_common()
        {
          tmr.print();
        }

	// ctor
	concurr_common(
	  const typename solver_t::rt_params_t &p,
          mem_t *mem_p,
	  const int &size
	) {
          // allocate the memory to be shared by multiple threads
          mem.reset(mem_p);
	  solver_t::alloc(mem.get(), p);

          // allocate per-thread structures
          init(p, p.grid_size, size); 
        }

        private:
 
        template <
          bcond::bcond_e type,
          bcond::drctn_e dir,
          int dim
        >
        void bcx_set(
          typename solver_t::bcp_t &bcp, 
          const std::array<int, solver_t::n_dims> &grid_size
        ) {
	  bcp.reset(
            new bcond::bcond<real_t, type, dir, solver_t::n_dims, dim>(
	      slab<dim>(grid_size), 
	      solver_t::halo, 
	      grid_size[0]
            )
          );
        }

        // 1D version
        void init(
          const typename solver_t::rt_params_t &p,
          const std::array<int, 1> &grid_size, const int &n0
        )
        {
          typename solver_t::bcp_t bxl, bxr, shrdl, shrdr;

          bcx_set<bcxl, bcond::left, 0>(bxl, grid_size);
	  bcx_set<bcxr, bcond::rght, 0>(bxr, grid_size);

	  for (int i0 = 0; i0 < n0; ++i0) 
          {
            shrdl.reset(new bcond::shared<real_t>());
            shrdr.reset(new bcond::shared<real_t>());

	    algos.push_back(
              new solver_t(
                typename solver_t::ctor_args_t({
		  mem.get(), 
		  i0 == 0      ? bxl : shrdl,
		  i0 == n0 - 1 ? bxr : shrdr,
		  slab<0>(grid_size, i0, n0)
                }), 
                p
              )
            );
          }
	}

        // 2D version
        // TODO: assert parallelisation in the right dimensions! (blitz::assertContiguous)
        void init(
          const typename solver_t::rt_params_t &p,
	  const std::array<int, 2> &grid_size, 
          const int &n0, const int &n1 = 1
        ) {
          for (int i0 = 0; i0 < n0; ++i0) 
          {
            for (int i1 = 0; i1 < n1; ++i1) 
            {
	      typename solver_t::bcp_t bxl, bxr, byl, byr, shrdl, shrdr;

              bcx_set<bcxl, bcond::left, 0>(bxl, grid_size);
	      bcx_set<bcxr, bcond::rght, 0>(bxr, grid_size);

              bcx_set<bcxl, bcond::left, 1>(byl, grid_size);
	      bcx_set<bcxr, bcond::rght, 1>(byr, grid_size);

              shrdl.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1
              shrdr.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1

              algos.push_back(
                new solver_t(
                  typename solver_t::ctor_args_t({
		    mem.get(), 
		    i0 == 0      ? bxl : shrdl,
		    i0 == n0 - 1 ? bxr : shrdr,
		    byl, byr, 
		    slab<0>(grid_size, i0, n0),  
                    slab<1>(grid_size, i1, n1)
                  }), 
                  p
                )
              );
            }
          }
	}

        // 3D version
        void init(
          const typename solver_t::rt_params_t &p,
	  const std::array<int, 3> &grid_size, 
          const int &n0, const int &n1 = 1, const int &n2 = 1
        ) {
          typename solver_t::bcp_t bxl, bxr, byl, byr, bzl, bzr, shrdl, shrdr;

	  // TODO: renew pointers only if invalid ?
	  for (int i0 = 0; i0 < n0; ++i0) 
          {
	    for (int i1 = 0; i1 < n1; ++i1) 
            {
	      for (int i2 = 0; i2 < n2; ++i2) 
              {
                bcx_set<bcxl, bcond::left, 0>(bxl, grid_size);
                bcx_set<bcxl, bcond::rght, 0>(bxr, grid_size);

                bcx_set<bcxl, bcond::left, 1>(byl, grid_size);
                bcx_set<bcxl, bcond::rght, 1>(byr, grid_size);

                bcx_set<bcxl, bcond::left, 2>(bzl, grid_size);
                bcx_set<bcxl, bcond::rght, 2>(bzr, grid_size);

                shrdl.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1
                shrdr.reset(new bcond::shared<real_t>()); // TODO: shrdy if n1 != 1

		algos.push_back(
                  new solver_t(
                    typename solver_t::ctor_args_t({
                      mem.get(), 
		      i0 == 0      ? bxl : shrdl,
		      i0 == n0 - 1 ? bxr : shrdr,
                      byl, byr, 
                      bzl, bzr, 
                      slab<0>(grid_size, i0, n0), 
                      slab<1>(grid_size, i1, n1), 
                      slab<2>(grid_size, i2, n2)
                    }), 
                    p
                  )
                );
              }
            }
          }
        }

        virtual void solve(int nt) = 0;

        public:
    
        void advance(int nt) final
        {   
          tmr.resume();
          solve(nt);
          tmr.stop();
        }  

	typename solver_t::arr_t advectee(int e = 0) final
	{
	  return mem->advectee(e);
	}

	typename solver_t::arr_t advector(int d = 0) final
	{
	  return mem->advector(d);
	}

	typename solver_t::arr_t g_factor() final
	{
	  return mem->g_factor();
	}

        bool *panic_ptr() final
        {
          return &this->mem->panic;
        }
      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace libmpdataxx
