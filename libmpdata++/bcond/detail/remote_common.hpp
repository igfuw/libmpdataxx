// common code for ``remote'' MPI boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

#if defined(USE_MPI)
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
      class remote_common : public detail::bcond_common<real_t, halo>
      {
	using parent_t = detail::bcond_common<real_t, halo>;

        protected:

	using arr_t = blitz::Array<real_t, n_dims>;
        using idx_t = blitz::RectDomain<n_dims>;

        private:

        const int grid_size_0;

#if defined(USE_MPI)
        boost::mpi::communicator mpicom;

#  if defined(NDEBUG)
        static const int n_reqs = 1; // data, reqs for recv only is enough?
#  else
        static const int n_reqs = 2; // data + ranges
#  endif

        std::array<boost::mpi::request, n_reqs> reqs;

        const int peer = dir == left
          ? (mpicom.rank() - 1 + mpicom.size()) % mpicom.size()
          : (mpicom.rank() + 1                ) % mpicom.size();
#endif

        protected:
        const bool is_cyclic = 
#if defined(USE_MPI)
	  (dir == left && mpicom.rank() == 0) ||
	  (dir == rght && mpicom.rank() == mpicom.size()-1);
#else
          false;
#endif

        void send(
          const arr_t &a, 
          const idx_t &idx_send 
        )
        {
          // distinguishing between left and right messages 
          // (important e.g. with 2 procs and cyclic bc)
          const int  
            msg_send = dir == left ? left : rght;

          // copying data to be sent (TODO: it doesn't work without copy(), why??)
	  arr_t buf_send = a(idx_send).copy();

#if defined(USE_MPI)
#  if !defined(NDEBUG)
          const int debug = 2;
#  endif
          // launching async data transfer
          {
            std::lock_guard<std::mutex> lock(libmpdataxx::concurr::detail::mpi_mutex);
	    mpicom.isend(peer, msg_send, buf_send);

            // sending debug information
#  if !defined(NDEBUG)
	    mpicom.isend(peer, msg_send ^ debug, std::pair<int,int>(
              idx_send[0].first(), 
              idx_send[0].last()
            ));
#  endif
#endif
          }
        };

        void recv(
          const arr_t &a, 
          const idx_t &idx_recv
        )
        {
          // distinguishing between left and right messages 
          // (important e.g. with 2 procs and cyclic bc)
          const int  
            msg_recv = dir == left ? rght : left;

          // copying data to be sent
	  arr_t 
            buf_recv(a(idx_recv).shape());

#if defined(USE_MPI)
#  if !defined(NDEBUG)
          const int debug = 2;
	  std::pair<int, int> buf_rng; 
#  endif
          // launching async data transfer
          {
            std::lock_guard<std::mutex> lock(libmpdataxx::concurr::detail::mpi_mutex);
	    reqs[0] = mpicom.irecv(peer, msg_recv, buf_recv);

            // sending debug information
#  if !defined(NDEBUG)
	    reqs[1] = mpicom.irecv(peer, msg_recv ^ debug, buf_rng);
#  endif
          }

          // waiting for the transfers to finish
	  boost::mpi::wait_all(reqs.begin(), reqs.end()); // MPI_Waitall is thread-safe?

          // checking debug information
          
          // positive modulo (grid_size_0 - 1)
          auto wrap = [this](int n) {return (n % (grid_size_0 - 1) + grid_size_0 - 1) % (grid_size_0 - 1);};

	  assert(wrap(buf_rng.first) == wrap(idx_recv[0].first()));
          assert(wrap(buf_rng.second) == wrap(idx_recv[0].last()));
#else
          assert(false);
#endif
          // writing received data to the array
	  a(idx_recv) = buf_recv;
        }


        void xchng(
          const arr_t &a, 
          const idx_t &idx_send, 
          const idx_t &idx_recv
        )
        {
          send(a, idx_send);
          recv(a, idx_recv);
        }

        public:

        // ctor                                  
        remote_common(                                                           
          const rng_t &i,
          const int &grid_size_0
        ) :
          parent_t(i, grid_size_0),
          grid_size_0(grid_size_0)
        {} 
      };
    }
  } // namespace bcond
} // namespace libmpdataxx
