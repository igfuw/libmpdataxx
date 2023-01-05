/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// memory management with (fractal) grid refinement

#pragma once

#include "sharedmem.hpp"

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
      class sharedmem_refined_common : public sharedmem<real_t, n_dims, n_tlev>
      {
        using parent_t = sharedmem<real_t, n_dims, n_tlev>;

        protected:

        using arr_t = typename parent_t::arr_t;

        blitz::TinyVector<int, n_dims> origin_ref;

        public:

        const int n_ref; // number of equal divisions of the large cell (in each direction), refined resolution is dx / n_ref;
                         // every n_ref scalar of the refined grid is at the same position as a scalar of the normal grid
                         // no refinement done in the halo, because there are no SD in the halo (it's not real space)
                         // what about MPI boundaries? there are refined points exactly at the boundary (since n_ref has to be even)
                         // note: if there is a refined cell that is divided by the MPI boundary, o we need to add contributions from microphysics to both processes on both sides?
                         //       maybe not, because microphysics contrbutions will affect the large cells, which are not divided by the MPI boundary...

        std::array<rng_t, n_dims> grid_size_ref;
        // TODO: these are public because used from outside in alloc - could friendship help?
        //arrvec_t<arr_t> GC_ref, psi_ref; 
        arrvec_t<arr_t> psi_ref; 

        // ctors
        sharedmem_refined_common(const std::array<int, n_dims> &grid_size, const int &size, const int &n_ref)
          : parent_t(grid_size, size), n_ref(n_ref)
        {
          assert(n_ref % 2 == 0 || n_ref == 1); // only division into even number of cells, because we assume that one of the refined scalar points is at the MPI boundary, which is in the middle between normal grid scalars

          if(n_ref % 2 == 0) // refinemenet actually done
          {
            // for now, require a grid_size that is convenient for fractal reconstruction (which calculates 2 points based on 3 points)
            // NOTE: fix this with proper halos (cyclic is easy, but what about rigid?)
            // NOTE2: this is actually a requirement for fractal reconstruction, not for any grid refinement, so move this somewhere else
            for (int d = 0; d < n_dims; ++d)
              if((grid_size[d] - 3) % 2 != 0) throw std::runtime_error("Fractal grid refinement requires nx/ny/nz = 3 + 2 * i, where i = 0,1,2,3,...");

            for (int d = 0; d < n_dims; ++d)
            {
              grid_size_ref[d] = refine_grid_size(
                this->grid_size[d],
                n_ref,
                d == 0 ? this->distmem.rank() : 0,
                d == 0 ? this->distmem.size() : 1
              );
              origin_ref[d] = grid_size_ref[d].first();

              this->distmem.grid_size_ref[d] = refine_grid_size(rng_t(0,grid_size[d]-1), n_ref, 0, 1).length();
            }
          }
          else if(n_ref == 1) // no refinement
          {
            for (int d = 0; d < n_dims; ++d)
            {
              grid_size_ref[d] = this->grid_size[d];
              origin_ref[d] = this->origin[d];
              this->distmem.grid_size_ref[d] = this->distmem.grid_size[d];
            }
          }
        }

        // NOTE: not all advectees are refined, so e (numbering) in refinee is different than in advectee
        virtual arr_t refinee(int e = 0) = 0;
      //  virtual const arr_t refinee_global_ref(int e = 0) = 0;

        public:
        static rng_t refine_grid_size(
          const rng_t &grid_size,
          const int &n_ref,
          const int &mpi_rank,
          const int &mpi_size
        ) {
          assert(n_ref % 2 == 0);
          // NOTE: overlapping points inbetween MPI domains
          return rng_t(
            mpi_rank == 0          ? grid_size.first() * n_ref : grid_size.first() * n_ref - n_ref / 2,
            mpi_rank == mpi_size-1 ? grid_size.last()  * n_ref : grid_size.last()  * n_ref + n_ref / 2 // refined points between MPI domains are evenly divided between MPI processes
          );
        }
      };



      template<typename real_t, int n_dims, int n_tlev>
      class sharedmem_refined
      {};

      template<typename real_t, int n_tlev>
      class sharedmem_refined<real_t, 3, n_tlev> : public sharedmem_refined_common<real_t, 3, n_tlev>
      {
        using parent_t = sharedmem_refined_common<real_t, 3, n_tlev>;
        using parent_t::parent_t; // inheriting ctors

        protected:
        using arr_t = typename parent_t::arr_t;

        public:
        arr_t refinee(int e = 0) override
        {
          return this->psi_ref[e](
            this->grid_size_ref[0],
            this->grid_size_ref[1],
            this->grid_size_ref[2]
          ).reindex(this->origin_ref);
        }
      };

    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
