/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

// ensuring thread-safe versions of system headers are used
#if !defined(_REENTRANT)
#  error _REENTRANT not defined, please use something like -pthread flag for gcc
#endif

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
#include <libmpdata++/bcond/polar_3d.hpp>
#include <libmpdata++/bcond/rigid_2d.hpp>
#include <libmpdata++/bcond/rigid_3d.hpp>
#include <libmpdata++/bcond/gndsky_3d.hpp>

#include <libmpdata++/concurr/detail/sharedmem.hpp>
#include <libmpdata++/concurr/detail/timer.hpp>
#include <libmpdata++/concurr/any.hpp>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      template<
        class solver_t_, 
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
      >
      class concurr_common : public any<typename solver_t_::real_t, solver_t_::n_dims, typename solver_t_::advance_arg_t>
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

	protected:

        // (cannot be nested due to templates)
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
        using advance_arg_t = typename solver_t::advance_arg_t;

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
	  solver_t::alloc(mem.get(), p.n_iters);

          // allocate per-thread structures
          init(p, mem->grid_size, size); 
        }

        private:
 
        template <
          bcond::bcond_e type,
          bcond::drctn_e dir,
          int dim
        >
        void bc_set(typename solver_t::bcp_t &bcp) 
        {
	  bcp.reset(
            new bcond::bcond<real_t, solver_t::halo, type, dir, solver_t::n_dims, dim>(
	      mem->slab(mem->grid_size[dim]), 
	      mem->grid_size[0].length() // TODO: get it from rt_params...
            )
          );
        }

        // 1D version
        void init(
          const typename solver_t::rt_params_t &p,
          const std::array<rng_t, 1> &grid_size, const int &n0
        )
        {
          typename solver_t::bcp_t bxl, bxr, shrdl, shrdr;

          bc_set<bcxl, bcond::left, 0>(bxl);
	  bc_set<bcxr, bcond::rght, 0>(bxr);

	  for (int i0 = 0; i0 < n0; ++i0) 
          {
            shrdl.reset(new bcond::shared<real_t, solver_t::halo>());
            shrdr.reset(new bcond::shared<real_t, solver_t::halo>());

	    algos.push_back(
              new solver_t(
                typename solver_t::ctor_args_t({
                  i0,
		  mem.get(), 
		  i0 == 0      ? bxl : shrdl,
		  i0 == n0 - 1 ? bxr : shrdr,
		  mem->slab(grid_size[0], i0, n0)
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
	  const std::array<rng_t, 2> &grid_size, 
          const int &n0, const int &n1 = 1
        ) {
          for (int i0 = 0; i0 < n0; ++i0) 
          {
            for (int i1 = 0; i1 < n1; ++i1) 
            {
	      typename solver_t::bcp_t bxl, bxr, byl, byr, shrdl, shrdr;

              bc_set<bcxl, bcond::left, 0>(bxl);
	      bc_set<bcxr, bcond::rght, 0>(bxr);

              bc_set<bcyl, bcond::left, 1>(byl);
	      bc_set<bcyr, bcond::rght, 1>(byr);

              shrdl.reset(new bcond::shared<real_t, solver_t::halo>()); // TODO: shrdy if n1 != 1
              shrdr.reset(new bcond::shared<real_t, solver_t::halo>()); // TODO: shrdy if n1 != 1

              algos.push_back(
                new solver_t(
                  typename solver_t::ctor_args_t({
                    i0,
		    mem.get(), 
		    i0 == 0      ? bxl : shrdl,
		    i0 == n0 - 1 ? bxr : shrdr,
		    byl, byr, 
		    mem->slab(grid_size[0], i0, n0),  
                    mem->slab(grid_size[1], i1, n1)
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
	  const std::array<rng_t, 3> &grid_size, 
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
                bc_set<bcxl, bcond::left, 0>(bxl);
                bc_set<bcxr, bcond::rght, 0>(bxr);

                bc_set<bcyl, bcond::left, 1>(byl);
                bc_set<bcyr, bcond::rght, 1>(byr);

                bc_set<bczl, bcond::left, 2>(bzl);
                bc_set<bczr, bcond::rght, 2>(bzr);

                shrdl.reset(new bcond::shared<real_t, solver_t::halo>()); // TODO: shrdy if n1 != 1
                shrdr.reset(new bcond::shared<real_t, solver_t::halo>()); // TODO: shrdy if n1 != 1

		algos.push_back(
                  new solver_t(
                    typename solver_t::ctor_args_t({
                      i0,
                      mem.get(), 
		      i0 == 0      ? bxl : shrdl,
		      i0 == n0 - 1 ? bxr : shrdr,
                      byl, byr, 
                      bzl, bzr, 
                      mem->slab(grid_size[0], i0, n0), 
                      mem->slab(grid_size[1], i1, n1), 
                      mem->slab(grid_size[2], i2, n2)
                    }), 
                    p
                  )
                );
              }
            }
          }
        }

        virtual void solve(advance_arg_t nt) = 0;

        public:
    
        void advance(advance_arg_t nt) final
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

        typename solver_t::arr_t vab_coefficient() final
	{
	  return mem->vab_coefficient();
	}
        
        typename solver_t::arr_t vab_relaxed_state(int d = 0) final
	{
	  return mem->vab_relaxed_state(d);
	}
	
        typename solver_t::arr_t sclr_array(const std::string &name, int n = 0) final
	{
	  return mem->sclr_array(name, n);
	}

        bool *panic_ptr() final
        {
          return &this->mem->panic;
        }

        const real_t time() const final
        {
          return algos[0].time_();
        }
      };
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
