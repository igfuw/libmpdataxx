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

#include <libmpdata++/concurr/detail/sharedmem.hpp>
#include <libmpdata++/concurr/detail/timer.hpp>
#include <libmpdata++/concurr/any.hpp>

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
#include <libmpdata++/bcond/remote_1d.hpp>
#include <libmpdata++/bcond/remote_2d.hpp>
#include <libmpdata++/bcond/remote_3d.hpp>
#include <libmpdata++/bcond/gndsky_3d.hpp>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      template <
        class real_t,
        class mem_t,
        bcond::drctn_e dir,
        int dim,
        int n_dims,
        int halo
      >
      void bc_set_remote(
        typename solver_t::bcp_t &bcp,
        const std::unique_ptr<mem_t> &mem;
        const int thread_rank,
        const int thread_count
      )
      {
        bcp.reset(
          new bcond::bcond<real_t, halo, bcond::remote, dir, n_dims, dim>(
            mem->slab(mem->grid_size[dim]),
            mem->distmem.grid_size
          )
        );
      }

      template <
        class real_t,
        class mem_t,
        bcond::drctn_e dir,
        int dim,
        int halo
      >
      void bc_set_remote<real_t, mem_t, dir, dim, 3, halo>(
        typename solver_t::bcp_t &bcp,
        const std::unique_ptr<mem_t> &mem;
        const int thread_rank,
        const int thread_count
      )
      {
        bcp.reset(
          new bcond::bcond<real_t, solver_t::halo, bcond::remote, dir, solver_t::n_dims, dim>(
            mem->slab(mem->grid_size[dim]),
            mem->distmem.grid_size,
            mem->slab(mem->grid_size[1], thread_rank, thread_count),
            thread_rank
          )
        );
      }

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
        void bc_set(
          typename solver_t::bcp_t &bcp,
          const int thread_rank  = 0, // required only by 3D remote (MPI) bcond
          const int thread_count = 0  // required only by 3D remote (MPI) bcond
        )
        {
          // sanity check - polar coords do not work with MPI yet
          if (type == bcond::polar && mem->distmem.size() > 1)
            throw std::runtime_error("Polar boundary conditions do not work with MPI.");

          // distmem overrides
          if (mem->distmem.size() > 1 && dim == 0)
          {
            if (
              // distmem domain interior
              (dir == bcond::left && mem->distmem.rank() > 0)
              ||
              (dir == bcond::rght && mem->distmem.rank() != mem->distmem.size() - 1)
              // cyclic condition for distmem domain (note: will not work if a non-cyclic condition is on the other end)
              ||
              (type == bcond::cyclic)
            )
            {
              // bc allocation, all mpi routines called by the remote bcnd ctor are thread-safe (?)
              bc_set_remote<real_t, mem_t, dir, dim, solver_t::n_dims, solver_t::halo>(
                new bcond::bcond<real_t, solver_t::halo, bcond::remote, dir, solver_t::n_dims, dim>(
                  mem->slab(mem->grid_size[dim]),
                  mem->distmem.grid_size,
                  mem->slab(mem->grid_size[1], thread_rank, thread_count),
                  thread_rank
                )
              );
              return;
            }
          }

          bcp.reset(
            new bcond::bcond<real_t, solver_t::halo, type, dir, solver_t::n_dims, dim>(
              mem->slab(mem->grid_size[dim]),
              mem->distmem.grid_size
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

          // NOTE: for remote bcond, thread_rank set to 0 on purpose in 1D to have propre left/right message tags 
          bc_set<bcxl, bcond::left, 0>(bxl);
          bc_set<bcxr, bcond::rght, 0>(bxr);

          for (int i0 = 0; i0 < n0; ++i0)
          {
            shrdl.reset(new bcond::shared<real_t, solver_t::halo, solver_t::n_dims>());
            shrdr.reset(new bcond::shared<real_t, solver_t::halo, solver_t::n_dims>());

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

              // NOTE: for remote bcond, thread_rank set to 0 on purpose in 2D to have propre left/right message tags 
              bc_set<bcxl, bcond::left, 0>(bxl);
              bc_set<bcxr, bcond::rght, 0>(bxr);

              bc_set<bcyl, bcond::left, 1>(byl);
              bc_set<bcyr, bcond::rght, 1>(byr);

              shrdl.reset(new bcond::shared<real_t, solver_t::halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1
              shrdr.reset(new bcond::shared<real_t, solver_t::halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1

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

        // 3D version, note sharedmem in y direction!
        void init(
          const typename solver_t::rt_params_t &p,
          const std::array<rng_t, 3> &grid_size,
          const int &n1, const int &n0 = 1, const int &n2 = 1
        ) {
          typename solver_t::bcp_t bxl, bxr, byl, byr, bzl, bzr, shrdl, shrdr;

          // TODO: renew pointers only if invalid ?
          for (int i0 = 0; i0 < n0; ++i0)
          {
            for (int i1 = 0; i1 < n1; ++i1)
            {
              for (int i2 = 0; i2 < n2; ++i2)
              {
                // i1 is the local thread rank, n1 is the number of threads. These are needed by remote bcond, because only rank=0 does mpi communication
                bc_set<bcxl, bcond::left, 0>(bxl, i1, n1);
                bc_set<bcxr, bcond::rght, 0>(bxr, i1, n1);

                bc_set<bcyl, bcond::left, 1>(byl);
                bc_set<bcyr, bcond::rght, 1>(byr);

                bc_set<bczl, bcond::left, 2>(bzl);
                bc_set<bczr, bcond::rght, 2>(bzr);

                shrdl.reset(new bcond::shared<real_t, solver_t::halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1
                shrdr.reset(new bcond::shared<real_t, solver_t::halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1

                algos.push_back(
                  new solver_t(
                    typename solver_t::ctor_args_t({
                      i1,
                      mem.get(),
                      bxl, bxr,
                      i1 == 0      ? byl : shrdl,
                      i1 == n1 - 1 ? byr : shrdr,
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

        const typename solver_t::arr_t advectee_global(int e = 0) final
        {
#if defined(USE_MPI)
          return mem->advectee_global(e);
#else
          return advectee(e);
#endif
        }

        void advectee_global_set(const typename solver_t::arr_t arr, int e = 0) final
        {
#if defined(USE_MPI)
          mem->advectee_global_set(arr, e);
#else
          advectee(e) = arr;
#endif
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

        const real_t min(int e = 0) const final
        {
          return mem->min(mem->advectee(e));
        }

        const real_t max(int e = 0) const final
        {
          return mem->max(mem->advectee(e));
        }
      };
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
