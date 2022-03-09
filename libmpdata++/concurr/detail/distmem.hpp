#pragma once

#if defined(USE_MPI)
#  include <boost/serialization/vector.hpp>
#  include <boost/mpi/communicator.hpp>
#  include <boost/mpi/collectives.hpp>
#else
#  include <cstdlib>
#endif

#include <numeric>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
#if defined(USE_MPI)
      namespace
      {
        bool mpi_initialized_before = true; // flag if MPI was initialized before distmem ctor

        // helpers that convert std::array into blitz::TinyVector, TODO: move them to formulas?
        blitz::TinyVector<int, 1> get_shape(std::array<int, 1> grid_size) {return blitz::shape(grid_size[0]);}
        blitz::TinyVector<int, 2> get_shape(std::array<int, 2> grid_size) {return blitz::shape(grid_size[0], grid_size[1]);}
        blitz::TinyVector<int, 3> get_shape(std::array<int, 3> grid_size) {return blitz::shape(grid_size[0], grid_size[1], grid_size[2]);}
      };
#endif

      template <typename real_t, int n_dims>
      class distmem
      {
#if defined(USE_MPI)
public: // TODO: just a temp measure, make it private again
        boost::mpi::communicator mpicom;
#endif

        private:

        template <typename Op, typename reduce_real_t> // some reductions done on different floating types (e.g. sum always on doubles)
        reduce_real_t reduce_hlpr(const reduce_real_t &val)
        {
#if defined(USE_MPI)
          reduce_real_t res;
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

        // sum always done on doubles
        // TODO: option to run more accurate summation (e.g. in pressure solver), see https://link.springer.com/content/pdf/10.1023/A:1008153532043.pdf
        double sum(const double &val)
        {
          return reduce_hlpr<std::plus<double>>(val);
        }

        template<class arr_t>
        const arr_t get_global_array(arr_t arr, const bool kij_to_kji)
        {
#if defined(USE_MPI)
          // a vector of number of elements to be sent by each non-root process
          std::vector<int> sizes(size());
          std::iota(sizes.begin(), sizes.end(), 0); // fill with 0,1,2,3,...

          for(auto &_size : sizes) 
          { 
            _size = domain_decomposition::slab(rng_t(0, grid_size[0]-1), _size, size()).length();
            for(int i=1; i<n_dims; ++i)
              _size *= grid_size[i];
          }

          // calc displacement
          std::vector<int> displ(sizes.size());
          std::partial_sum(sizes.begin(), sizes.end(), displ.begin());
          std::transform(displ.begin(), displ.end(), sizes.begin(), displ.begin(), std::minus<int>()); // exclusive_scan is c++17
          // a vector that will store the received data, relevant only on process rank=0
          std::vector<real_t> out_values(std::accumulate(grid_size.begin(), grid_size.end(), 1, std::multiplies<int>()));

          std::vector<real_t> in_values_vec;

          if(kij_to_kji)
          {
            // create an array that will store advectee to be sent in a contiguous memory block using the (default) kji storage order
            // NOTE: libmpdata++ 3d blitz arrays, like advectee, are in the kij order
            arr_t in_values_arr(arr.shape());
            in_values_arr = arr;
            // wrap in_values_arr in a std::vector
            in_values_vec = std::vector<real_t>(in_values_arr.begin(), in_values_arr.end());
          }
          else
          {
            // create an array that will store arr to be sent in a contiguous memory block
            in_values_vec = std::vector<real_t>(arr.size());
            std::copy(arr.begin(), arr.end(), in_values_vec.begin());
          }

          // gather the data from all processes on rank=0
          boost::mpi::gatherv(mpicom, in_values_vec, out_values.data(), sizes, displ, 0);
          // send the result to other processes
          boost::mpi::broadcast(mpicom, out_values, 0);

          blitz::Array<real_t, n_dims> res(out_values.data(), get_shape(grid_size), blitz::duplicateData);
          return res;
#else
          return arr;
#endif
        }

        // ctor
        distmem(const std::array<int, n_dims> &grid_size)
          : grid_size(grid_size)
        {
#if !defined(USE_MPI)
          if (
            // mvapich2
            std::getenv("MV2_COMM_WORLD_RANK") != NULL && std::atoi(std::getenv("MV2_COMM_WORLD_SIZE")) > 1 ||
            // mpich
            std::getenv("PMI_RANK") != NULL && std::atoi(std::getenv("PMI_SIZE")) > 1 ||
            // openmpi
            std::getenv("OMPI_COMM_WORLD_RANK") != NULL && std::atoi(std::getenv("OMPI_COMM_WORLD_SIZE")) > 1 ||
            // lam
            std::getenv("LAMRANK") != NULL && std::atoi(std::getenv("LAMSIZE")) > 1
          ) throw std::runtime_error("mpirun environment variable detected but libmpdata++ was compiled with MPI disabled");
#else
          // init mpi here, since distmem is constructed before hdf5
          // will be finalized in slvr_common dtor, since distmem is destructed before hdf5;
          // boost::mpi::enviro being global var caused deadlocks when freeing memory
          int mpi_initialized;
          MPI_Initialized(&mpi_initialized);
          if(!mpi_initialized)
          {
            mpi_initialized_before = false;
            int th_lvl_provided;
            MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &th_lvl_provided);
          }
          if (boost::mpi::environment::thread_level() != boost::mpi::threading::multiple)
          {
            throw std::runtime_error("failed to initialise MPI environment with MPI_THREAD_MULTIPLE");
          }
          mpicom = boost::mpi::communicator(MPI_COMM_WORLD, boost::mpi::comm_duplicate); // use a duplicate of MPI_COMM_WORLD, can't construct it before MPI_Init call (?)
#endif
        }
      };
    }
  }
}
