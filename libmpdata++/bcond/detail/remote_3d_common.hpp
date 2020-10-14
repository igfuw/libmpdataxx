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

        protected:

        using arr_t = typename parent_t::arr_t;
        using idx_t = typename parent_t::idx_t;

        const int thread_rank, thread_size;

        private:

        const rng_t thread_j;
        const int grid_size_y;

        // try to guess what should be the whole domain exchanged by this process
        // based on the difference between idx to be sent by this thread and idx of this process
        idx_t extend_idx(idx_t idx)
        {
//std::cerr << "extend idx start idx(1): " << idx.lbound(1) << ", " << idx.ubound(1) << std::endl;
          idx.lbound(1) = 0                 + idx.lbound(1) - thread_j.first(); // does it have to start at 0?
          idx.ubound(1) = (grid_size_y - 1) + idx.ubound(1) - thread_j.last();  // does it have to end at grid_size_y - 1?
//std::cerr << "extend idx end idx(1): " << idx.lbound(1) << ", " << idx.ubound(1) << std::endl;
          return idx;
        }

        public:

        void xchng (
          const arr_t &a,
          const idx_t &idx_send,
          const idx_t &idx_recv
        )
        {
          parent_t::xchng(a, extend_idx(idx_send), extend_idx(idx_recv));
        }

        void send (
          const arr_t &a,
          const idx_t &idx_send
        )
        {
          parent_t::send(a, extend_idx(idx_send));
        }

        void recv (
          const arr_t &a,
          const idx_t &idx_recv
        )
        {
          parent_t::recv(a, extend_idx(idx_recv));
        }

        // ctor
        remote_3d_common(
          const rng_t &i,
          const std::array<int, 3> &distmem_grid_size,
          const rng_t _thread_j,
          const int thread_rank,
          const int thread_size
        ) :
          parent_t(i, distmem_grid_size),
          thread_rank(thread_rank),
          thread_size(thread_size),
          thread_j(_thread_j),
          grid_size_y(distmem_grid_size[1])
        {
#if defined(USE_MPI)
          // only 2 threads do mpi, others don't need buffers
          if(thread_rank != 0 && thread_rank != thread_size-1)
          {
            free(parent_t::buf_send);
            free(parent_t::buf_recv);
          }
//std::cerr << "remote_3d_common ctor thread_j: " << thread_j.lbound(0) << ", " << thread_j.ubound(0) << std::endl;
//std::cerr << "remote_3d_common ctor _thread_j: " << _thread_j.lbound(0) << ", " << _thread_j.ubound(0) << std::endl;
//std::cerr << "remote_3d_common ctor f-l thread_j: " << thread_j.first() << ", " << thread_j.last() << std::endl;
//std::cerr << "remote_3d_common ctor f-l _thread_j: " << _thread_j.first() << ", " << _thread_j.last() << std::endl;
#endif
        }

        // dtor
        ~remote_3d_common()
        {
#if defined(USE_MPI)
          if(thread_rank == 0 || thread_rank == thread_size-1)
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
