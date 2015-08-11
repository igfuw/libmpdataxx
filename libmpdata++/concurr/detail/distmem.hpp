#pragma once

#if defined(USE_MPI)
#  include <boost/mpi/environment.hpp>
#  include <boost/mpi/communicator.hpp>
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
        boost::mpi::environment mpienv;
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
      };
    }
  }
}
