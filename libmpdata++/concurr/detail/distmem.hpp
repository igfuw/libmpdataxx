#pragma once

#if defined(USE_MPI)
#  include <boost/mpi/communicator.hpp>
#  include <boost/mpi/collectives.hpp>
#else
#  include <cstdlib>
#endif


namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
#if defined(USE_MPI)
      bool mpi_initialized_before = true; // flag if MPI was initialized before distmem ctor
#endif

      template <typename real_t, int n_dims>
      class distmem
      {
#if defined(USE_MPI)
public: // TODO: just a temp measure, make it private again
        boost::mpi::communicator mpicom;
#endif

        private:

        template <typename Op>
        real_t reduce_hlpr(const real_t &val)
        {
#if defined(USE_MPI)
          real_t res;
          boost::mpi::all_reduce(mpicom, val, res, Op()); // it's thread-safe?
          return res;
#else
          return val;
#endif
        }


        public:

        std::array<int, n_dims> grid_size;

        int rank() 
        { 
#if defined(USE_MPI)
          return mpicom.rank();  // is it thread-safe? TODO: init once in ctor?
#else
          return 0;
#endif
        }

        int size() 
        {
#if defined(USE_MPI)
          return mpicom.size();   // is it thread-safe? TODO: init once in ctor?
#else
          return 1;
#endif
        }

        void barrier()
        {
#if defined(USE_MPI)
          mpicom.barrier();
#else
          assert(false);
#endif
        }

        real_t min(const real_t &val)
        {
#if defined(USE_MPI)
          return reduce_hlpr<boost::mpi::minimum<real_t>>(val);
#else
          return val;
#endif
        }

        real_t max(const real_t &val)
        {
#if defined(USE_MPI)
          return reduce_hlpr<boost::mpi::maximum<real_t>>(val);
#else
          return val;
#endif
        }

        real_t sum(const real_t &val)
        {
          return reduce_hlpr<std::plus<real_t>>(val);
        }

        // ctor
        distmem(const std::array<int, n_dims> &grid_size) 
          : grid_size(grid_size)
        {
#if !defined(USE_MPI)
          if (
            // mvapich2
            std::getenv("MV2_COMM_WORLD_RANK") != NULL ||
            // mpich
            std::getenv("PMI_RANK") != NULL ||
            // openmpi
            std::getenv("OMPI_COMM_WORLD_RANK") != NULL ||
            // lam
            std::getenv("LAMRANK") != NULL
          ) throw std::runtime_error("mpirun environment variable detected but libmpdata++ was compiled with MPI disabled");
#else
          // init mpi here, since distmem is constructed before hdf5
          // will be finalized in slvr_common dtor, since distmem is destructed before hdf5;
          // boost::mpi::enviro being global var caused deadlocks when freeing memory
          if(!MPI::Is_initialized())
          {
            mpi_initialized_before = false;
            MPI::Init_thread(MPI_THREAD_MULTIPLE);
          }
          if (boost::mpi::environment::thread_level() != boost::mpi::threading::multiple)
          {
            throw std::runtime_error("failed to initialise MPI environment with MPI_THREAD_MULTIPLE");
          }
#endif
        }
      };
    }
  }
}
