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
      class distmem
      {
#if defined(USE_MPI)
        // note: this is intentionally left commented out
        //       as it is not used in libmpdata++, it seems unneeded, and it causes:
        //       "*** The MPI_Errhandler_set() function was called after MPI_FINALIZE was invoked."
        //       error if multiple solvers are instantiated 
        // boost::mpi::environment mpienv; 

        boost::mpi::communicator mpicom;
#endif

        public:

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

        // ctor
        distmem() 
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
