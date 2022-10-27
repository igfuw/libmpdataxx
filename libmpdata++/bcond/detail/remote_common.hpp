// common code for ``remote'' MPI boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

#if defined(USE_MPI)
#  include <boost/serialization/vector.hpp>
#  include <boost/mpi/communicator.hpp>
#  include <boost/mpi/nonblocking.hpp>
#endif

namespace libmpdataxx
{
  namespace bcond
  {
    namespace detail
    {
      template <typename real_t, int halo, drctn_e dir, int n_dims>
      class remote_common : public detail::bcond_common<real_t, halo, n_dims>
      {
        using parent_t = detail::bcond_common<real_t, halo, n_dims>;

        protected:

        using arr_t = blitz::Array<real_t, n_dims>;
        using idx_t = blitz::RectDomain<n_dims>;

        real_t *buf_send,
               *buf_recv;

        private:

#if defined(USE_MPI)
        boost::mpi::communicator mpicom;

#  if defined(NDEBUG)
        static const int n_reqs = 2; // data, reqs for recv only is enough?
        static const int n_dbg_send_reqs = 0;
        static const int n_dbg_tags = 0;
#  else
        static const int n_reqs = 4; // data + ranges
        static const int n_dbg_send_reqs = 1;
        static const int n_dbg_tags = 2;
#  endif

        std::array<boost::mpi::request, n_reqs> reqs;

        const int peer = dir == left
          ? (mpicom.rank() - 1 + mpicom.size()) % mpicom.size()
          : (mpicom.rank() + 1                ) % mpicom.size();

#  if !defined(NDEBUG)
          std::pair<int, int> buf_rng;
#  endif
#endif

        protected:
        const bool is_cyclic =
#if defined(USE_MPI)
          (dir == left && mpicom.rank() == 0) ||
          (dir == rght && mpicom.rank() == mpicom.size()-1);
#else
          false;
#endif

        void send_hlpr(
          const arr_t &a,
          const idx_t &idx_send
        )
        {
#if defined(USE_MPI)
          // distinguishing between left and right messages
          // (important e.g. with 2 procs and cyclic bc)
          const int
            msg_send = dir == left ? left : rght;

          std::cerr << "send_hlpr idx dir " << dir << " : " 
            << " (" << idx_send.lbound(0) << ", " << idx_send.ubound(0) << ")"  
            << " (" << idx_send.lbound(1) << ", " << idx_send.ubound(1) << ")"  
            << " (" << idx_send.lbound(2) << ", " << idx_send.ubound(2) << ")"  
            << std::endl;

          std::cerr << "buf_send: " << buf_send << std::endl;

          // arr_send references part of the send buffer that will be used
          arr_t arr_send(buf_send, a(idx_send).shape(), blitz::neverDeleteData);
          // copying data to be sent
          arr_send = a(idx_send);
          std::cerr << "arr_send: " << arr_send << std::endl;

          // launching async data transfer
          if(arr_send.size()!=0)
          {
            // use the pointer+size kind of send instead of serialization of blitz arrays, because
            // serialization caused memory leaks, probably because it breaks blitz reference counting
            reqs[0] = mpicom.isend(peer, msg_send, buf_send, arr_send.size());

            // sending debug information
#  if !defined(NDEBUG)
            reqs[1] = mpicom.isend(peer, msg_send + n_dbg_tags, std::pair<int,int>(
              idx_send[0].first(),
              idx_send[0].last()
            ));
#  endif
          }
#else
          assert(false);
#endif
        };

        void recv_hlpr(
          const arr_t &a,
          const idx_t &idx_recv
        )
        {
#if defined(USE_MPI)
          const int
            msg_recv = dir == left ? rght : left;

          std::cerr << "recv_hlpr idx dir " << dir << " : " 
            << " (" << idx_recv.lbound(0) << ", " << idx_recv.ubound(0) << ")"  
            << " (" << idx_recv.lbound(1) << ", " << idx_recv.ubound(1) << ")"  
            << " (" << idx_recv.lbound(2) << ", " << idx_recv.ubound(2) << ")"  
            << std::endl;


          // launching async data transfer
          if(a(idx_recv).size()!=0) // TODO: test directly size of idx_recv
          {
            reqs[1+n_dbg_send_reqs] = mpicom.irecv(peer, msg_recv, buf_recv, a(idx_recv).size());

            // sending debug information
#  if !defined(NDEBUG)
            reqs[3] = mpicom.irecv(peer, msg_recv + n_dbg_tags, buf_rng);
#  endif
          }
#else
          assert(false);
#endif
        }

        void send(
          const arr_t &a,
          const idx_t &idx_send
        )
        {
#if defined(USE_MPI)
          send_hlpr(a, idx_send);

          // waiting for the transfers to finish
          boost::mpi::wait_all(reqs.begin(), reqs.begin() + 1 + n_dbg_send_reqs); // MPI_Waitall is thread-safe?
#else
          assert(false);
#endif
        }

        void recv(
          const arr_t &a,
          const idx_t &idx_recv
        )
        {
#if defined(USE_MPI)
          //auto arr_recv = recv_hlpr(a, idx_recv);
          recv_hlpr(a, idx_recv);

          // waiting for the transfers to finish
          boost::mpi::wait_all(reqs.begin() + 1 + n_dbg_send_reqs, reqs.end()); // MPI_Waitall is thread-safe?

          // a blitz handler for the used part of the receive buffer
          arr_t arr_recv(buf_recv, a(idx_recv).shape(), blitz::neverDeleteData); // TODO: shape directly from idx_recv

          // checking debug information

          // positive modulo (grid_size_0 - 1)
//          auto wrap = [this](int n) {return (n % (grid_size_0 - 1) + grid_size_0 - 1) % (grid_size_0 - 1);};
//          assert(wrap(buf_rng.first) == wrap(idx_recv[0].first()));
//          assert(wrap(buf_rng.second) == wrap(idx_recv[0].last()));

          // writing received data to the array
          a(idx_recv) = arr_recv;
#else
          assert(false);
#endif
        }

        void xchng(
          const arr_t &a,
          const idx_t &idx_send,
          const idx_t &idx_recv
        )
        {
#if defined(USE_MPI)
          send_hlpr(a, idx_send);
          recv_hlpr(a, idx_recv);

          // waiting for the transfers to finish
          boost::mpi::wait_all(reqs.begin(), reqs.end());

          // a blitz handler for the used part of the receive buffer
          arr_t arr_recv(buf_recv, a(idx_recv).shape(), blitz::neverDeleteData);

          // checking debug information

          // positive modulo (grid_size_0 - 1)
         // auto wrap = [this](int n) {return (n % (grid_size_0 - 1) + grid_size_0 - 1) % (grid_size_0 - 1);};
         // assert(wrap(buf_rng.first) == wrap(idx_recv[0].first()));
         // assert(wrap(buf_rng.second) == wrap(idx_recv[0].last()));

          // writing received data to the array
          a(idx_recv) = arr_recv;
#else
          assert(false);
#endif
        }

        public:

        // ctor
        remote_common(
          boost::mpi::communicator &mpic,
          const rng_t &i,
          const std::array<int, n_dims> &distmem_grid_size,
          bool single_threaded = false,
          const int thread_rank = -1, 
          const int thread_size = -1
        ) :
          parent_t(i, distmem_grid_size, single_threaded, thread_rank, thread_size)
          //,mpicom(mpic, boost::mpi::comm_take_ownership)
          ,mpicom(mpic, boost::mpi::comm_attach)
        {
#if defined(USE_MPI)
  
          const int slice_size = n_dims==1 ? 1 : (n_dims==2? distmem_grid_size[1]+6 : (distmem_grid_size[1]+6) * (distmem_grid_size[2]+6) ); // 3 is the max halo size (?), so 6 on both sides
std::cerr << "remote_common ctor, " 
  << " distmem_grid_size[0]: " << distmem_grid_size[0]
  << " distmem_grid_size[1]: " << distmem_grid_size[1]
  << " distmem_grid_size[2]: " << distmem_grid_size[2]
  << " slice_size: " << slice_size 
  << " halo: " << halo
  << std::endl;
          // allocate enough memory in buffers to store largest halos to be sent
          buf_send = (real_t *) malloc(halo * slice_size * sizeof(real_t));
          buf_recv = (real_t *) malloc(halo * slice_size * sizeof(real_t));
#endif
        }

        // dtor
        ~remote_common()
        {
#if defined(USE_MPI)
          free(buf_send);
          free(buf_recv);
#endif
        }
      };
    }
  } // namespace bcond
} // namespace libmpdataxx

