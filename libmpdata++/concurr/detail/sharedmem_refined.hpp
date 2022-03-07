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

        const int n_ref; // number of equal divisions of the large cell (in each direction)

        std::array<rng_t, n_dims> grid_size_ref;
        // TODO: these are public because used from outside in alloc - could friendship help?
        arrvec_t<arr_t> GC_ref, psi_ref; 

        // ctors
        sharedmem_refined_common(const std::array<int, n_dims> &grid_size, const int &size, const int &n_ref)
          : parent_t(grid_size, size), n_ref(n_ref)
        {
          for (int d = 0; d < n_dims; ++d)
          {
            grid_size_ref[d] = refine_grid_size(this->grid_size[d], n_ref);
            origin_ref[d] = grid_size_ref[d].first();
            this->distmem.grid_size_ref[d] = grid_size_ref[d].last() - grid_size_ref[d].first() + 1;
          }
        }

        virtual arr_t refinee(int e = 0) = 0;
      //  virtual const arr_t refinee_global_ref(int e = 0) = 0;

        public:
        static rng_t refine_grid_size(
          const rng_t &grid_size,
          const int &n_ref
        ) {
          return rng_t(
            grid_size.first() * n_ref,
            (grid_size.last() + 1) * n_ref - 1
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
          std::cerr << "refinee: " << this->psi_ref[e](
            this->grid_size_ref[0],
            this->grid_size_ref[1],
            this->grid_size_ref[2]
          ).reindex(this->origin_ref);

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
