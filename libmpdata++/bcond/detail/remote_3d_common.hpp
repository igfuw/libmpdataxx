// common code for ``remote'' MPI boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/remote_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    namespace detail
    {
      template <typename real_t, int halo, drctn_e dir>
      class remote_3d_common : public remote_common<real_t, halo, dir, 3>
      {
        using parent_t = detail::remote_common<real_t, halo, dir, 3>;

        using arr_t = blitz::Array<real_t, 3>;
        using idx_t = blitz::RectDomain<3>;

        const int thread_rank;
        const rng_t thread_j;
        const int grid_size_y;

        // try to guess what should be the whole domain exchanged by this process
        // based on the difference between idx to be sent by this thread and idx of this process
        idx_t extend_idx(idx_t idx)
        {
          idx.lbound(1) = 0 + idx.lbound(1) - thread_j.lbound(0); // does it have to start at 0?
          idx.ubound(1) = grid_size_y + idx.ubound(1) - thread_j.ubound(0); // does it have to end at grid_size_y? what about courants?
          return idx;
        }

        public:

        void xchng (
          const arr_t &a,
          const idx_t &idx_send,
          const idx_t &idx_recv
        ) override
        {
          if(thread_rank != 0) return;
          parent_t::xchng(a, extend_idx(idx_send), extend_idx(idx_recv));
        }

        void send (
          const arr_t &a,
          const idx_t &idx_send
        ) override
        {
          if(thread_rank != 0) return;
          parent_t::send(a, extend_idx(idx_send));
        }

        void recv (
          const arr_t &a,
          const idx_t &idx_recv
        ) override
        {
          if(thread_rank != 0) return;
          parent_t::recv(a, extend_idx(idx_recv));
        }

        // ctor
        remote_3d_common(
          const rng_t &i,
          const std::array<int, 3> &distmem_grid_size,
          const rng_t thread_j,
          const int thread_rank
        ) :
          parent_t(i, distmem_grid_size),
          thread_rank(thread_rank),
          thread_j(thread_j),
          grid_size_y(distmem_grid_size[1])
        {
#if defined(USE_MPI)
          // only thread 0 does mpi, others don't need buffers
          if(thread_rank != 0)
          {
            free(parent_t::buf_send);
            free(parent_t::buf_recv);
          }
#endif
        }

        // dtor
        ~remote_3d_common()
        {
#if defined(USE_MPI)
          if(thread_rank == 0)
          {
            free(parent_t::buf_send);
            free(parent_t::buf_recv);
          }
          // hack to make free in ~remote_common give defined behaviour
          parent_t::buf_send = nullptr;
          parent_t::buf_recv = nullptr;
#endif
        }
      };
    };
  } // namespace bcond
} // namespace libmpdataxx
