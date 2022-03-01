/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// memory management with (fractal) grid refinement

#pragma once

#include <sharedmem.hpp>

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
      class sharedmem_refined_common :: sharedmem<real_t, n_dims, n_tlev>
      {
        using parent_t = sharedmem<real_t, n_dims, n_tlev>;
        using arr_t = typename parent_t::parent_t::arr_t;

        protected:

        const int n_ref; // number of equal divisions of the large cell (in each direction)
        blitz::TinyVector<int, n_dims> origin_ref;

        public:

        std::array<rng_t, n_dims> grid_size_ref;
        // TODO: these are public because used from outside in alloc - could friendship help?
        arrvec_t<arr_t> GC_ref, psi_ref; 

        // ctors
        sharedmem_refined_common(const std::array<int, n_dims> &grid_size, const int &size, const int &n_ref)
          : parent_t(grid_size, size), n_ref(n_ref)
        {
          for (int d = 0; d < n_dims; ++d)
          {
            grid_size_ref[d] = refine_grid_size(this->grid_size, n_ref)
            origin_ref[d] = grid_size_ref[d].first();
          }
        }

        virtual arr_t advectee_ref(int e = 0) = 0;
        virtual const arr_t advectee_global_ref(int e = 0) = 0;

        public:
        static rng_t refine_grid_size(
          const rng_t &grid_size,
          const int &n_ref
        ) {
          return rng_t(
            grid_size.first() * n_ref,
            grid_size.first() + (grid_size.last() - grid_size.first()) * n_ref
          );
        }
      };








/*




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
          arr_t *ret = new arr_t(arg->dataFirst(), arg->shape(), blitz::neverDeleteData, blitz::GeneralArrayStorage<3>(arg->ordering(), {true, true, true}));
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
            // create an array that will store advectee to be sent in a contiguous memory block using the (default) kji storage order
            // NOTE: libmpdata++ 3d blitz arrays, like advectee, are in the kij order
            blitz::Array<real_t, 3> in_values_arr(advectee(e).shape());
            in_values_arr = advectee(e);
            // wrap in_values_arr in a std::vector
            std::vector<real_t> in_values_vec(in_values_arr.begin(), in_values_arr.end());

            // gather the data from all processes on rank=0
            boost::mpi::gatherv(this->distmem.mpicom, in_values_vec, out_values.data(), sizes, displ, 0);
            // send the result to other processes
            boost::mpi::broadcast(this->distmem.mpicom, out_values, 0);

            blitz::Array<real_t, 3> res(out_values.data(), blitz::shape(
              this->distmem.grid_size[0], this->grid_size[1].length(), this->grid_size[2].length()),
              blitz::duplicateData
              );
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
      */
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
