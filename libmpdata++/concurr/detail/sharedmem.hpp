/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <unordered_map>
#include <boost/ptr_container/ptr_vector.hpp>

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/formulae/arakawa_c.hpp>
#include <libmpdata++/concurr/detail/distmem.hpp>

#include <array>
#include <numeric>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      template <
        typename real_t,
        int n_dims,
        int n_tlev
      >
      class sharedmem_common
      {

        static_assert(n_dims > 0, "n_dims <= 0");
        static_assert(n_tlev > 0, "n_tlev <= 0");

        std::unique_ptr<blitz::Array<real_t, 1>> xtmtmp;
        std::unique_ptr<blitz::Array<double, 1>> sumtmp;

        protected:

        using arr_t = blitz::Array<real_t, n_dims>;
        blitz::TinyVector<int, n_dims> origin;

        public:

        int n = 0;
        const int size;
        std::array<rng_t, n_dims> grid_size;
        bool panic = false; // for multi-threaded SIGTERM handling

        // dimension in which sharedmem domain decomposition is done
        // 1D and 2D - domain decomposed in 0-th dimension (x)
        // 3D - domain decomposed in 1-st dimension (y) for better workload balance in MPI runs (MPI is decomposed in x)
        const int shmem_decomp_dim;

        detail::distmem<real_t, n_dims> distmem;

        // TODO: these are public because used from outside in alloc - could friendship help?
        arrvec_t<arr_t> GC, ndt_GC, ndtt_GC;
        std::vector<arrvec_t<arr_t>> psi; // TODO: since n_eqns is known, could make it an std::array!
        std::unique_ptr<arr_t> G;
        std::unique_ptr<arr_t> vab_coeff; // velocity absorber coefficient
        arrvec_t<arr_t> vab_relax; // velocity absorber relaxed state
        arrvec_t<arr_t> khn_tmp; // Kahan sum for donor-cell

        std::unordered_map<
          const char*, // intended for addressing with __FILE__
          boost::ptr_vector<arrvec_t<arr_t>>
        > tmp;

        // list of temporary fields that can be accessed from outside of concurr
        std::unordered_map<
          std::string,
          std::pair<const char*, int>
        > avail_tmp;

        virtual void barrier()
        {
          assert(false && "sharedmem_common::barrier() called!");
        }

        void cycle(const int &rank)
        {
          barrier();
          if (rank == 0) n = (n + 1) % n_tlev - n_tlev; // - n_tlev assures Python-type end-of-array cyclic behaviour works
          barrier();
        }

        // ctors
        // TODO: fill reducetmp with NaNs (or use 1-element arrvec_t - it's NaN-filled by default)
        sharedmem_common(const std::array<int, n_dims> &grid_size, const int &size)
          : n(0), distmem(grid_size), size(size), shmem_decomp_dim(n_dims < 3 ? 0 : 1) // TODO: is n(0) needed?
        {
          for (int d = 0; d < n_dims; ++d)
          {
            this->grid_size[d] = slab(
              rng_t(0, grid_size[d]-1),
              d == 0 ? distmem.rank() : 0,          // decomposition along x, because that's MPI decomposition
              d == 0 ? distmem.size() : 1
             // d == shmem_decomp_dim ? distmem.rank() : 0,
             // d == shmem_decomp_dim ? distmem.size() : 1
            );
            origin[d] = this->grid_size[d].first();
          }

          if (size > grid_size[0])
            throw std::runtime_error("number of subdomains greater than number of gridpoints");

          if (n_dims != 1)
            sumtmp.reset(new blitz::Array<double, 1>(this->grid_size[n_dims-2])); 
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }

        /// @brief concurrency-aware summation of array elements
        double sum(const int &rank, const arr_t &arr, const idx_t<n_dims> &ijk, const bool sum_khn)
        {
          // doing a two-step sum to reduce numerical error
          // and make parallel results reproducible
          for (int c = ijk[shmem_decomp_dim].first(); c <= ijk[shmem_decomp_dim].last(); ++c) // TODO: optimise for i.count() == 1
          {
            auto slice_idx = ijk;
            slice_idx.lbound(shmem_decomp_dim) = c;
            slice_idx.ubound(shmem_decomp_dim) = c;

            if (sum_khn)
              (*sumtmp)(c) = blitz::kahan_sum(arr(slice_idx));
            else
              (*sumtmp)(c) = blitz::sum(arr(slice_idx));
          }
          barrier(); // wait for all threads to calc their part
#if !defined(USE_MPI)
          double result;
          if (sum_khn)
            result = blitz::kahan_sum(*sumtmp);
          else
            result = blitz::sum(*sumtmp);
          barrier();
          return result;
#else
          if(rank == 0)
          {
            // master thread calculates the sum from this process, stores in shared array
            if (sum_khn)
              (*sumtmp)(grid_size[shmem_decomp_dim].first())= blitz::kahan_sum(*sumtmp); // inplace?!
            else
              (*sumtmp)(grid_size[shmem_decomp_dim].first())= blitz::sum(*sumtmp); // inplace?!
            // master thread calculates sum of sums from all processes
            (*sumtmp)(grid_size[shmem_decomp_dim].first()) = this->distmem.sum((*sumtmp)(grid_size[shmem_decomp_dim].first())); // inplace?!
          }
          barrier();
          double res = (*sumtmp)(grid_size[shmem_decomp_dim].first()); // propagate the total sum to all threads of the process
          barrier(); // to avoid sumtmp being overwritten by next call to sum from other thread
          return res;
#endif
        }

        /// @brief concurrency-aware summation of a (element-wise) product of two arrays
        double sum(const int &rank, const arr_t &arr1, const arr_t &arr2, const idx_t<n_dims> &ijk, const bool sum_khn)
        {
          // doing a two-step sum to reduce numerical error
          // and make parallel results reproducible
          for (int c = ijk[shmem_decomp_dim].first(); c <= ijk[shmem_decomp_dim].last(); ++c)
          {
            auto slice_idx = ijk;
            slice_idx.lbound(shmem_decomp_dim) = c;
            slice_idx.ubound(shmem_decomp_dim) = c;

            if (sum_khn)
              (*sumtmp)(c) = blitz::kahan_sum(arr1(slice_idx) * arr2(slice_idx));
            else
              (*sumtmp)(c) = blitz::sum(arr1(slice_idx) * arr2(slice_idx));
          }
          // TODO: code below same as in the function above
          barrier(); // wait for all threads to calc their part
#if !defined(USE_MPI)
          double result;
          if (sum_khn)
            result = blitz::kahan_sum(*sumtmp);
          else
            result = blitz::sum(*sumtmp);
          barrier();
          return result;
#else
          if(rank == 0)
          {
            // master thread calculates the sum from this process, stores in shared array
            if (sum_khn)
              (*sumtmp)(grid_size[shmem_decomp_dim].first())= blitz::kahan_sum(*sumtmp); // inplace?!
            else
              (*sumtmp)(grid_size[shmem_decomp_dim].first())= blitz::sum(*sumtmp); // inplace?!
            // master thread calculates sum of sums from all processes
            (*sumtmp)(grid_size[shmem_decomp_dim].first()) = this->distmem.sum((*sumtmp)(grid_size[shmem_decomp_dim].first())); // inplace?!
          }
          barrier();
          double res = (*sumtmp)(grid_size[shmem_decomp_dim].first()); // propagate the total sum to all threads of the process
          barrier(); // to avoid sumtmp being overwritten by next call to sum from other thread
          return res;
#endif
        }

        real_t min(const int &rank, const arr_t &arr)
        {
          // min across local threads
          (*xtmtmp)(rank) = blitz::min(arr);
          barrier();
#if !defined(USE_MPI)
          real_t result = blitz::min(*xtmtmp);
          barrier();
          return result;
#else
          if(rank == 0)
          {
            (*xtmtmp)(0) = blitz::min(*xtmtmp);
            // min across mpi processes
            (*xtmtmp)(0) = this->distmem.min((*xtmtmp)(0));
          }
          barrier();
          real_t res = (*xtmtmp)(0); // propagate the total min to all threads of the process
          barrier(); // to avoid xtmtmp being overwritten by some other threads' next sum call
          return res;
#endif
        }

        // TODO: almost the same as min
        real_t max(const int &rank, const arr_t &arr)
        {
          // max across local threads
          (*xtmtmp)(rank) = blitz::max(arr);
          barrier();
#if !defined(USE_MPI)
          real_t result = blitz::max(*xtmtmp);
          barrier();
          return result;
#else
          if(rank == 0)
          {
            (*xtmtmp)(0) = blitz::max(*xtmtmp);
            // max across mpi processes
            (*xtmtmp)(0) = this->distmem.max((*xtmtmp)(0));
          }
          barrier();
          real_t res = (*xtmtmp)(0); // propagate the total max to all threads of the process
          barrier(); // to avoid xtmtmp being overwritten by some other threads' next sum call
          return res;
#endif
        }

        // single-threaded, MPI-aware versions of the min and max functions
        real_t min(const arr_t &arr)
        {
          // min across local threads
          real_t result = blitz::min(arr);
          // min across mpi processes
          result = this->distmem.min(result);
          return result;
        }

        real_t max(const arr_t &arr)
        {
          // min across local threads
          real_t result = blitz::max(arr);
          // min across mpi processes
          result = this->distmem.max(result);
          return result;
        }

        // this hack is introduced to allow to use neverDeleteData
        // and hence to not use BZ_THREADSAFE
        private:
        boost::ptr_vector<arr_t> tobefreed;

        public:
        virtual arr_t *never_delete(arr_t *arg)
        {
          arr_t *ret = new arr_t(arg->dataFirst(), arg->shape(), blitz::neverDeleteData);
          ret->reindexSelf(arg->base());
          return ret;
        }

        arr_t *old(arr_t *arg)
        {
          tobefreed.push_back(arg);
          arr_t *ret = this->never_delete(arg);
          return ret;
        }

        private:
        // helper methods to define subdomain ranges
        static int min(const int &span, const int &rank, const int &size)
        {
          return rank * span / size;
        }

        static int max(const int &span, const int &rank, const int &size)
        {
          return min(span, rank + 1, size) - 1;
        }

        public:
        static rng_t slab(
          const rng_t &span,
          const int &rank = 0,
          const int &size = 1
        ) {
          return rng_t(
            span.first() + min(span.length(), rank, size),
            span.first() + max(span.length(), rank, size)
          );
        }

        virtual arr_t advectee(int e = 0) = 0;

        void advectee_global_set(const arr_t arr, int e = 0)
        {
#if defined(USE_MPI)
          if(this->distmem.size() > 1)
          {
            advectee(e) = arr(slab(rng_t(0, distmem.grid_size[0]-1), distmem.rank(), distmem.size()));
          }
          else
#endif
          advectee(e) = arr;
        }

        protected:

        rng_t distmem_ext(const rng_t &rng)
        {
          return rng_t(rng.first()-1, rng.last());
        }

      };

      template<typename real_t, int n_dims, int n_tlev>
      class sharedmem
      {};

      template<typename real_t, int n_tlev>
      class sharedmem<real_t, 1, n_tlev> : public sharedmem_common<real_t, 1, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 1, n_tlev>;
        using parent_t::parent_t; // inheriting ctors

        public:

        // accessor methods
        blitz::Array<real_t, 1> advectee(int e = 0)
        {
          assert(this->n < n_tlev);

          // returning just the domain interior, i.e. without halos
          // reindexing so that element 0 is at 0
          return this->psi[e][ this->n ](
            this->grid_size[0]
          ).reindex(this->origin);
        }

        const blitz::Array<real_t, 1> advectee_global(int e = 0)
        {
#if defined(USE_MPI)
          if(this->distmem.size() > 1)
          {
// TODO: move some of that to distmem...
// TODO: a lot of common code betwee 1,2 and 3 dimensions...

            // a vector of number of elements to be sent by each non-root process
            std::vector<int> sizes(this->distmem.size());
            std::iota(sizes.begin(), sizes.end(), 0); // fill with 1,2,3,...
            for(auto &size : sizes) { size = this->slab(rng_t(0, this->distmem.grid_size[0]-1), size, this->distmem.size()).length();}
            // calc displacement
            std::vector<int> displ(sizes.size());
            std::partial_sum(sizes.begin(), sizes.end(), displ.begin());
            std::transform(displ.begin(), displ.end(), sizes.begin(), displ.begin(), std::minus<int>()); // exclusive_scan is c++17
            // a vector that will store the received data, relevant only on process rank=0
            std::vector<real_t> out_values(this->distmem.grid_size[0]);
            // create an array that will store advectee to be sent in a contiguous memory block
            std::vector<real_t> in_values_vec(advectee(e).size());
            std::copy(advectee(e).begin(), advectee(e).end(), in_values_vec.begin());

            // gather the data from all processes on rank=0
            boost::mpi::gatherv(this->distmem.mpicom, in_values_vec, out_values.data(), sizes, displ, 0);
            // send the result to other processes
            boost::mpi::broadcast(this->distmem.mpicom, out_values, 0);

            blitz::Array<real_t, 1> res(out_values.data(), blitz::shape(this->distmem.grid_size[0]), blitz::duplicateData);
            return res;
          }
          else
#endif
            return advectee(e);
        }

        blitz::Array<real_t, 1> advector(int d = 0)
        {
          using namespace arakawa_c;
          assert(d == 0);
          // returning just the domain interior, i.e. without halos
          // reindexed to make it more intuitive when working with index placeholders
          // (i.e. border between cell 0 and cell 1 is indexed with 0)
          auto orgn = decltype(this->origin)({
                 this->origin[0] - 1
               });

          return this->GC[d](
            this->distmem_ext(this->grid_size[0]^(-1)^h)
          ).reindex(
            this->distmem.rank() > 0
              ? decltype(this->origin)({this->origin[0] - 1})
              : orgn
          );
        }

        blitz::Array<real_t, 1> g_factor()
        {
          // a sanity check
          if (this->G.get() == nullptr)
            throw std::runtime_error("g_factor() called with nug option unset?");

          // the same logic as in advectee() - see above
          return (*this->G)(
            this->grid_size[0]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 1> vab_coefficient()
        {
          throw std::logic_error("absorber not yet implemented in 1d");
        }

        blitz::Array<real_t, 1> vab_relaxed_state(int d = 0)
        {
          throw std::logic_error("absorber not yet implemented in 1d");
        }

        blitz::Array<real_t, 1> sclr_array(const std::string& name, int n = 0)
        {
          return this->tmp.at(this->avail_tmp[name].first)[this->avail_tmp[name].second][n](
            this->grid_size[0]
          ).reindex(this->origin);
        }
      };

      template<typename real_t, int n_tlev>
      class sharedmem<real_t, 2, n_tlev> : public sharedmem_common<real_t, 2, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 2, n_tlev>;
        using parent_t::parent_t; // inheriting ctors

        public:

        blitz::Array<real_t, 2> advectee(int e = 0)
        {
          assert(this->n < n_tlev);

          return this->psi[e][ this->n ](
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
        }

        const blitz::Array<real_t, 2> advectee_global(int e = 0)
        {
#if defined(USE_MPI)
          if(this->distmem.size() > 1)
          {
// TODO: move some of that to distmem...
// TODO: a lot of common code betwee 1,2 and 3 dimensions...

            // a vector of number of elements to be sent by each non-root process
            std::vector<int> sizes(this->distmem.size());
            std::iota(sizes.begin(), sizes.end(), 0); // fill with 1,2,3,...
            for(auto &size : sizes)
            {
              size = this->slab(rng_t(0, this->distmem.grid_size[0]-1), size, this->distmem.size()).length()
                      * this->grid_size[1].length();
            }
            // calc displacement
            std::vector<int> displ(sizes.size());
            std::partial_sum(sizes.begin(), sizes.end(), displ.begin());
            std::transform(displ.begin(), displ.end(), sizes.begin(), displ.begin(), std::minus<int>()); // exclusive_scan is c++17
            // a vector that will store the received data, relevant only on process rank=0
            std::vector<real_t> out_values(this->distmem.grid_size[0] * this->grid_size[1].length());
            // create an array that will store advectee to be sent in a contiguous memory block
            std::vector<real_t> in_values_vec(advectee(e).size());
            std::copy(advectee(e).begin(), advectee(e).end(), in_values_vec.begin());

            // gather the data from all processes on rank=0
            boost::mpi::gatherv(this->distmem.mpicom, in_values_vec, out_values.data(), sizes, displ, 0);
            // send the result to other processes
            boost::mpi::broadcast(this->distmem.mpicom, out_values, 0);

            blitz::Array<real_t, 2> res(out_values.data(), blitz::shape(
              this->distmem.grid_size[0], this->grid_size[1].length()),
              blitz::duplicateData);
            return res;
          }
          else
#endif
            return advectee(e);
        }

        blitz::Array<real_t, 2> advector(int d = 0)
        {
          using namespace arakawa_c;
          assert(d == 0 || d== 1);
          // returning just the domain interior, i.e. without halos
          // reindexed to make it more intuitive when working with index placeholders
          auto orgn = decltype(this->origin)({
                 this->origin[0] - 1,
                 this->origin[1]
               });

          switch (d)
          {
            case 0:
              return this->GC[d](
                this->distmem_ext(this->grid_size[0]^(-1)^h),
                this->grid_size[1]
              ).reindex(orgn);
            case 1:
              return this->GC[d](
                this->distmem_ext(this->grid_size[0]),
                this->grid_size[1]^(-1)^h
              ).reindex(orgn);
            default: assert(false); throw;
          }
        }

        blitz::Array<real_t, 2> g_factor()
        {
          // a sanity check
          if (this->G.get() == nullptr)
            throw std::runtime_error("g_factor() called with nug option unset?");

          // the same logic as in advectee() - see above
          return (*this->G)(
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 2> vab_coefficient()
        {
          // a sanity check
          if (this->vab_coeff.get() == nullptr)
            throw std::runtime_error("vab_coeff() called with option vip_vab unset?");

          // the same logic as in advectee() - see above
          return (*this->vab_coeff)(
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 2> vab_relaxed_state(int d = 0)
        {
          assert(d == 0 || d== 1);
          // a sanity check
          if (this->vab_coeff.get() == nullptr)
            throw std::runtime_error("vab_relaxed_state() called with option vip_vab unset?");
          // the same logic as in advectee() - see above
          return this->vab_relax[d](
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 2> sclr_array(const std::string& name, int n = 0)
        {
          return this->tmp.at(this->avail_tmp[name].first)[this->avail_tmp[name].second][n](
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
        }
      };

      template<typename real_t, int n_tlev>
      class sharedmem<real_t, 3, n_tlev> : public sharedmem_common<real_t, 3, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 3, n_tlev>;
        using arr_t = typename parent_t::arr_t;
        using parent_t::parent_t; // inheriting ctors

        public:

        virtual arr_t *never_delete(arr_t *arg) override
        {
          arr_t *ret = new arr_t(arg->dataFirst(), arg->shape(), blitz::neverDeleteData, arr3D_storage);
          ret->reindexSelf(arg->base());
          return ret;
        }

        blitz::Array<real_t, 3> advectee(int e = 0)
        {
          assert(this->n < n_tlev);

          return this->psi[e][ this->n ](
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
        }

        const blitz::Array<real_t, 3> advectee_global(int e = 0)
        {
#if defined(USE_MPI)
          if(this->distmem.size() > 1)
          {
// TODO: move some of that to distmem...
// TODO: a lot of common code betwee 1,2 and 3 dimensions...

            // a vector of number of elements to be sent by each non-root process
            std::vector<int> sizes(this->distmem.size());
            std::iota(sizes.begin(), sizes.end(), 0); // fill with 1,2,3,...
            for(auto &size : sizes)
            {
              size = this->slab(rng_t(0, this->distmem.grid_size[0]-1), size, this->distmem.size()).length()
                      * this->grid_size[1].length() * this->grid_size[2].length();
            }
            // calc displacement
            std::vector<int> displ(sizes.size());
            std::partial_sum(sizes.begin(), sizes.end(), displ.begin());
            std::transform(displ.begin(), displ.end(), sizes.begin(), displ.begin(), std::minus<int>()); // exclusive_scan is c++17
            // a vector that will store the received data, relevant only on process rank=0
            std::vector<real_t> out_values(this->distmem.grid_size[0] * this->grid_size[1].length() * this->grid_size[2].length());
            // create an array that will store advectee to be sent in a contiguous memory block
            std::vector<real_t> in_values_vec(advectee(e).size());
            std::copy(advectee(e).begin(), advectee(e).end(), in_values_vec.begin());

            // gather the data from all processes on rank=0
            boost::mpi::gatherv(this->distmem.mpicom, in_values_vec, out_values.data(), sizes, displ, 0);
            // send the result to other processes
            boost::mpi::broadcast(this->distmem.mpicom, out_values, 0);

            blitz::Array<real_t, 3> res(out_values.data(), blitz::shape(
              this->distmem.grid_size[0], this->grid_size[1].length(), this->grid_size[2].length()),
              blitz::duplicateData);
            return res;
          }
          else
#endif
            return advectee(e);
        }

        blitz::Array<real_t, 3> advector(int d = 0)
        {
          using namespace arakawa_c;
          assert(d == 0 || d == 1 || d == 2);
          // returning just the domain interior, i.e. without halos
          // reindexed to make it more intuitive when working with index placeholders
          auto orgn = decltype(this->origin)({
                 this->origin[0] - 1,
                 this->origin[1],
                 this->origin[2]
               });

          switch (d)
          {
            case 0:
              return this->GC[d](
                this->distmem_ext(this->grid_size[0]^(-1)^h),
                this->grid_size[1],
                this->grid_size[2]
              ).reindex(orgn);
            case 1:
              return this->GC[d](
                this->distmem_ext(this->grid_size[0]),
                this->grid_size[1]^(-1)^h,
                this->grid_size[2]
              ).reindex(orgn);
            case 2:
              return this->GC[d](
                this->distmem_ext(this->grid_size[0]),
                this->grid_size[1],
                this->grid_size[2]^(-1)^h
              ).reindex(orgn);
            default: assert(false); throw;
          }
        }

        blitz::Array<real_t, 3> g_factor()
        {
          // a sanity check
          if (this->G.get() == nullptr)
            throw std::runtime_error("g_factor() called with nug option unset?");

          // the same logic as in advectee() - see above
          return (*this->G)(
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 3> vab_coefficient()
        {
          // a sanity check
          if (this->vab_coeff.get() == nullptr)
            throw std::runtime_error("vab_coeff() called with option vip_vab unset?");

          // the same logic as in advectee() - see above
          return (*this->vab_coeff)(
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 3> vab_relaxed_state(int d = 0)
        {
          assert(d == 0 || d == 1 || d == 2);
          // a sanity check
          if (this->vab_coeff.get() == nullptr)
            throw std::runtime_error("vab_relaxed_state() called with option vip_vab unset?");
          // the same logic as in advectee() - see above
          return this->vab_relax[d](
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 3> sclr_array(const std::string& name, int n = 0)
        {
          return this->tmp.at(this->avail_tmp[name].first)[this->avail_tmp[name].second][n](
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
        }
      };
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
