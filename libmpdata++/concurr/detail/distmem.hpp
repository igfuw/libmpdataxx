#pragma once

#if defined(USE_MPI)
#  include <boost/mpi/environment.hpp>
#  include <boost/mpi/communicator.hpp>
#else
#  include <cstdlib>
#endif


namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      namespace 
      {
#if defined(USE_MPI)
        // the anonymous namespace is a hack to avoid
        // "*** The MPI_Errhandler_set() function was called after MPI_FINALIZE was invoked."
        // error if multiple solvers are instantiated 
        boost::mpi::environment mpienv;
#endif
      }

      template <int n_dims>
      class distmem
      {
#if defined(USE_MPI)
        boost::mpi::communicator mpicom;
#endif

        public:

        std::array<int, n_dims> grid_size;

        int rank() 
        { 
#if defined(USE_MPI)
          return mpicom.rank(); 
#else
          return 0;
#endif
        }

        int size() 
        {
#if defined(USE_MPI)
          return mpicom.size(); 
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

        // ctor
        distmem(const std::array<int, n_dims> &grid_size) 
          : grid_size(grid_size)
        {
#if !defined(USE_MPI)
          if (
            // mpich
            std::getenv("PMI_RANK") != NULL ||
            // openmpi
            std::getenv("OMPI_COMM_WORLD_RANK") != NULL ||
            // lam
            std::getenv("LAMRANK") != NULL
          ) throw std::runtime_error("mpirun environment variable detected but libmpdata++ was compiled with MPI disabled");
#endif
        }
      };
    }
  }
}
