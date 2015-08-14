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
      template <typename real_t, int halo, drctn_e dir>
      class remote_common : public detail::bcond_common<real_t, halo>
      {
	using parent_t = detail::bcond_common<real_t, halo>;
	using arr_t = blitz::Array<real_t, 1>;

        const int grid_size_0;

        // zero-base Blitz++ arrays used as buffers
        arr_t 
          buf_send = arr_t(halo), 
          buf_recv = arr_t(halo);

#if defined(USE_MPI)
        boost::mpi::communicator mpicom;

#  if defined(NDEBUG)
        static const int n_reqs = 2; // data
#  else
        static const int n_reqs = 4; // data + ranges
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

        void xchng(
          const arr_t &a, 
          const rng_t &rng_send, 
          const rng_t &rng_recv
        )
        {
          // checking if the buffers are long enough
          assert(rng_send.length() == buf_send.extent(0));
          assert(rng_recv.length() == buf_recv.extent(0));

          // distinguishing between left and right messages 
          // (important e.g. with 2 procs and cyclic bc)
          const int  
            msg_send = dir == left ? left : rght,
            msg_recv = dir == left ? rght : left;

          // copying data to be sent 
          // (couldn't sort out how to manage non-zero-base Blitz++ arrays with Boost.MPI)
	  buf_send = a(rng_send);

#if defined(USE_MPI)
          // launching async data transfer
	  reqs[0] = mpicom.isend(peer, msg_send, buf_send);
	  reqs[1] = mpicom.irecv(peer, msg_recv, buf_recv);

          // sending debug information
#  if !defined(NDEBUG)
          const int debug = 2;
	  reqs[2] = mpicom.isend(peer, msg_send ^ debug, std::pair<int,int>(
            rng_send.first(), 
            rng_send.last()
          ));
	  std::pair<int, int> buf_rng; 
	  reqs[3] = mpicom.irecv(peer, msg_recv ^ debug, buf_rng);
#  endif

          // waiting for the transfers to finish
	  boost::mpi::wait_all(reqs.begin(), reqs.end());

          // checking debug information
	  assert(buf_rng.first  % grid_size_0 == rng_recv.first() % grid_size_0);
          assert(buf_rng.second % grid_size_0 == rng_recv.last()  % grid_size_0);
#else
          assert(false);
#endif
          // writing received data to the array
	  a(rng_recv) = buf_recv;
        }

        public:

        // ctor                                  
        remote_common(                                                           
          const rng_t &i,
          const int grid_size_0
        ) :
          parent_t(i, grid_size_0),
          grid_size_0(grid_size_0)
        {} 
      };
    }
  } // namespace bcond
} // namespace libmpdataxx
