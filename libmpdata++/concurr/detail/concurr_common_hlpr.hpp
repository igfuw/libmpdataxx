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

#include <libmpdata++/concurr/detail/sharedmem_refined.hpp>
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
#if defined(USE_MPI)
  #include <libmpdata++/bcond/remote_1d.hpp>
  #include <libmpdata++/bcond/remote_2d.hpp>
  #include <libmpdata++/bcond/remote_3d.hpp>
#endif
#include <libmpdata++/bcond/gndsky_3d.hpp>

#include <libmpdata++/solvers/detail/solver_type_traits.hpp>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      // helpers for setting remote bcond
#if defined(USE_MPI)
      template <
        class real_t,
        bcond::drctn_e dir,
        int dim,
        int n_dims,
        int halo,
        class bcp_t,
        class mem_t,
        class grid_size_t,
        class distmem_grid_size_t
      >
      struct bc_set_remote_impl
      {
        static void _(
          bcp_t &bcp,
          const std::unique_ptr<mem_t> &mem,
          const grid_size_t grid_size,
          const distmem_grid_size_t distmem_grid_size,
          const int thread_rank,
          const int thread_size,
          boost::mpi::communicator &mpic
        )
        {
          bcp.reset(
            new bcond::bcond<real_t, halo, bcond::remote, dir, n_dims, dim>(
              mpic,
              mem->slab(grid_size[dim]),
              distmem_grid_size,
              thread_rank,
              thread_size
            )
          );
        }
      };

      // 3d specialization
      template <
        class real_t,
        bcond::drctn_e dir,
        int dim,
        int halo,
        class bcp_t,
        class mem_t,
        class grid_size_t,
        class distmem_grid_size_t
      >
      struct bc_set_remote_impl<real_t, dir, dim, 3, halo, bcp_t, mem_t, grid_size_t, distmem_grid_size_t>
      {
        static void _(
          bcp_t &bcp,
          const std::unique_ptr<mem_t> &mem,
          const grid_size_t grid_size,
          const distmem_grid_size_t distmem_grid_size,
          const int thread_rank,
          const int thread_size,
          boost::mpi::communicator &mpic
        )
        {
          bcp.reset(
            new bcond::bcond<real_t, halo, bcond::remote, dir, 3, dim>(
              mpic,
              mem->slab(grid_size[dim]),
              distmem_grid_size,
              mem->slab(grid_size[1], thread_rank, thread_size), // NOTE: we assume here remote 3d bcond is only on the edges perpendicular to x
              thread_rank,
              thread_size
            )
          );
        }
      };

      template <
        class real_t,
        bcond::drctn_e dir,
        int dim,
        int n_dims,
        int halo,
        class bcp_t,
        class mem_t,
        class grid_size_t,
        class distmem_grid_size_t
      >
      void bc_set_remote(
        bcp_t &bcp,
        const std::unique_ptr<mem_t> &mem,
        const grid_size_t grid_size,
        const distmem_grid_size_t distmem_grid_size,
        const int thread_rank,
        const int thread_size,
        boost::mpi::communicator &mpic
      )
      {
        bc_set_remote_impl<real_t, dir, dim, n_dims, halo, bcp_t, mem_t, grid_size_t, distmem_grid_size_t>::_(bcp, mem, grid_size, distmem_grid_size, thread_rank, thread_size, mpic);
      }
#endif

      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
      >
      class concurr_common_hlpr : public any<typename solver_t_::real_t, solver_t_::n_dims, typename solver_t_::advance_arg_t>
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

        using mem_t = typename solver_t::mem_t;

        // member fields
        boost::ptr_vector<solver_t> algos;
        std::unique_ptr<mem_t> mem;
        timer tmr;

        typename solver_t::bcp_t bxl, bxr, byl, byr, bzl, bzr, shrdl, shrdr;

        public:

        typedef typename solver_t::real_t real_t;
        using advance_arg_t = typename solver_t::advance_arg_t;

        // dtor
        virtual ~concurr_common_hlpr()
        {
          tmr.print();
        }

        // ctor
        concurr_common_hlpr(
          const typename solver_t::rt_params_t &p,
          mem_t *mem_p,
          const int &size
        ) {
          // allocate the memory to be shared by multiple threads
          mem.reset(mem_p);
          solver_t::alloc(mem.get(), p.n_iters);
        }

        protected:

        template <
          bcond::bcond_e type,
          bcond::drctn_e dir,
          int dim,
          int halo,
          class bcp_t,
          class grid_size_t,
          class distmem_grid_size_t,
          class mpicomm_t
        >
        void bc_set(
          bcp_t &bcp,
          const grid_size_t grid_size,
          const distmem_grid_size_t distmem_grid_size,
          mpicomm_t &mpic,
          const int thread_rank  = 0,  // required only by 3D remote (MPI) and open bconds
          const int thread_size = 0    // required only by 3D remote (MPI) and open bconds
        )
        {
          // sanity check - polar coords do not work with MPI yet
          if (type == bcond::polar && mem->distmem.size() > 1)
            throw std::runtime_error("Polar boundary conditions do not work with MPI.");

          // distmem overrides
#if defined(USE_MPI)
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
              bc_set_remote<real_t, dir, dim, solver_t::n_dims, halo>(
                bcp,
                mem,
                grid_size,
                distmem_grid_size,
                thread_rank,
                thread_size,
                mpic
              );
              return;
            }
          }
#endif

          // 3d open bcond needs to know thread rank and size, because it zeroes perpendicular vectors
          if (type == bcond::open && solver_t::n_dims == 3)
          {
            bcp.reset(
              new bcond::bcond<real_t, halo, type, dir, solver_t::n_dims, dim>(
                mem->slab(grid_size[dim]),
                distmem_grid_size,
                false,
                thread_rank,
                thread_size
              )
            );
            return;
          }

          // else: not remote and not open_3d
          bcp.reset(
            new bcond::bcond<real_t, halo, type, dir, solver_t::n_dims, dim>(
              mem->slab(grid_size[dim]),
              distmem_grid_size
            )
          );
        }

#if defined(USE_MPI)
        template<int halo, class bcp_t, class grid_size_t, class distmem_grid_size_t>
        void init_bcs_3d(
          bcp_t &bxl, bcp_t &bxr, bcp_t &byl, bcp_t &byr, bcp_t &bzl, bcp_t &bzr, bcp_t &shrdl, bcp_t &shrdr,
          const grid_size_t &grid_size,
          const distmem_grid_size_t &distmem_grid_size,
          const int &i1, const int &n1,
          boost::mpi::communicator &mpic)
        {
          // i1 is the local thread rank, n1 is the number of threads. These are needed by remote bcond, because only rank=0 does mpi communication
          //boost::mpi::communicator mpic(MPI_COMM_WORLD, boost::mpi::comm_duplicate);
          //boost::mpi::communicator mpic(MPI_COMM_WORLD, boost::mpi::comm_attach);

          bc_set<bcxl, bcond::left, 0, halo>(bxl, grid_size, distmem_grid_size, mpic, i1, n1);
          bc_set<bcxr, bcond::rght, 0, halo>(bxr, grid_size, distmem_grid_size, mpic, i1, n1);
                                                                                    
          bc_set<bcyl, bcond::left, 1, halo>(byl, grid_size, distmem_grid_size, mpic);
          bc_set<bcyr, bcond::rght, 1, halo>(byr, grid_size, distmem_grid_size, mpic);
                                                                                    
          bc_set<bczl, bcond::left, 2, halo>(bzl, grid_size, distmem_grid_size, mpic);
          bc_set<bczr, bcond::rght, 2, halo>(bzr, grid_size, distmem_grid_size, mpic);

          shrdl.reset(new bcond::shared<real_t, halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1
          shrdr.reset(new bcond::shared<real_t, halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1
        }

#else

        template<int halo, class bcp_t, class grid_size_t, class distmem_grid_size_t>
        void init_bcs_3d(
          bcp_t &bxl, bcp_t &bxr, bcp_t &byl, bcp_t &byr, bcp_t &bzl, bcp_t &bzr, bcp_t &shrdl, bcp_t &shrdr,
          const grid_size_t &grid_size,
          const distmem_grid_size_t &distmem_grid_size,
          const int &i1, const int &n1)
        {
          // i1 is the local thread rank, n1 is the number of threads. These are needed by remote bcond, because only rank=0 does mpi communication
          bc_set<bcxl, bcond::left, 0, halo>(bxl, grid_size, distmem_grid_size, i1, n1);
          bc_set<bcxr, bcond::rght, 0, halo>(bxr, grid_size, distmem_grid_size, i1, n1);

          bc_set<bcyl, bcond::left, 1, halo>(byl, grid_size, distmem_grid_size, i1, n1);
          bc_set<bcyr, bcond::rght, 1, halo>(byr, grid_size, distmem_grid_size, i1, n1);

          bc_set<bczl, bcond::left, 2, halo>(bzl, grid_size, distmem_grid_size, i1, n1);
          bc_set<bczr, bcond::rght, 2, halo>(bzr, grid_size, distmem_grid_size, i1, n1);

          shrdl.reset(new bcond::shared<real_t, halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1
          shrdr.reset(new bcond::shared<real_t, halo, solver_t::n_dims>()); // TODO: shrdy if n1 != 1
        }
#endif

        private:
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

      template< class mem_t, class solver_t>
      mem_t* mem_factory(const typename solver_t::rt_params_t &p)
      {
        if constexpr (solvers::detail::slvr_with_frac_recn<typename solver_t::ct_params_t_>())
          return new mem_t(p.grid_size, pow(2, p.n_fra_iter));
        else
          return new mem_t(p.grid_size);
      }
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
