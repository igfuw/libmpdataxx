/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

// solver with fractal reconstruction of SGS fields

#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>
//#include <numeric>
//#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t, int minhalo>
      class mpdata_rhs_vip_prs_sgs_fra_common : public mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
      {
        using parent_t = mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>;

        public:

        using real_t = typename ct_params_t::real_t;

        protected:

//        const int n_fra; // number of fields with fractal reconstruction
        const int n_ref, // number of refinements; refined resolution is dx / n_ref
                  n_fra_iter;
        
        // refined arrays have no halos (do they need them? halos can be added by extending grid_size in alloc() in mpdata_rhs_vip_prs_sgs_fra.hpp by e.g. mem->n_ref/2)
//        const int halo_ref; // size of halos of refined scalar (and vector?) arrays

        // TODO: make these const!
        idx_t<ct_params_t::n_dims>  ijk_ref; // range of refinee handled by given solver
        const idxs_t<ct_params_t::n_dims> ijk_r2r; // resolved to refined; refined scalars at the same position as resolved scalars

//        static rng_t rng_sclr_ref(const rng_t &rng) { return rng^halo_ref; }

        static rng_t rng_midpoints(const rng_t &rng, const int rank = 0, const int size = 1, const bool overlap = true) // input range, rank and size of thread / MPI process, should created ranges have overlaps between different ranks 
        {
          assert(rng.stride() % 2 == 0);
          if(rng.last() - rng.first() < rng.stride()) return rng;
          else if (overlap)
            return rng_t(
              rank == 0 ? rng.first() + rng.stride() / 2 : rng.first() - rng.stride() / 2,
              rank == size - 1 ? rng.last() - rng.stride() / 2 : rng.last() + rng.stride() / 2, 
              rng.stride()); 
          else
            return rng_t(
              rng.first() + rng.stride() / 2,
              rank == size - 1 ? rng.last() - rng.stride() / 2 : rng.last() + rng.stride() / 2, 
              rng.stride()); 
        }

        static rng_t rng_midpoints_out(const rng_t &rng, const int rank = 0, const int size = 1) 
        {
          assert(rng.stride() % 4 == 0);
          if(rng.last() - rng.first() < rng.stride()) return rng;
          else return rng_t(
            rank == 0 ? rng.first() - rng.stride() / 4 : rng.first() + rng.stride() / 4,
            rank == size - 1 ? rng.last() + rng.stride() / 4 : rng.last() - rng.stride() / 4, 
            rng.stride() / 2); 
        }

        static rng_t rng_half_stride(const rng_t &rng, const int rank = 0, const int size = 1, const bool overlap = true) 
        {
          assert(rng.stride() % 2 == 0);
          if(overlap)
            return rng_t(
              rank > 0 ? rng.first() - rng.stride() / 2 : rng.first(), 
              rank < size - 1 ? rng.last() + rng.stride() / 2 : rng.last(), 
              rng.stride() / 2); 
          else
            return rng_t(
              rng.first(), 
              rank < size - 1 ? rng.last() + rng.stride() / 2 : rng.last(), 
              rng.stride() / 2); 
        }

        static rng_t rng_dbl_stride(const rng_t &rng) 
        {
          return rng_t(rng.first(), rng.last(), 2*rng.stride());
        }

        public:

        struct rt_params_t : parent_t::rt_params_t
        {
          unsigned int n_fra_iter = 1; // number of iterations of fractal reconstruction
        };

        // ctor
        mpdata_rhs_vip_prs_sgs_fra_common(
          typename parent_t::ctor_args_t args,
          const rt_params_t &p
        ) :
          parent_t(args, p),
          n_fra_iter(p.n_fra_iter),
          n_ref(this->mem->n_ref),
//          halo_ref(n_ref / 2),
          // ijk_ref init below assumes 3D (shmem decomp dim along y);
          // TODO: move this to concurr_common::init()? add something to ctor_args_t?
          ijk_r2r{
            {this->ijk[0].first() * n_ref, this->ijk[1].first() * n_ref, this->ijk[2].first() * n_ref}, // lbound
            {this->ijk[0].last() * n_ref, this->ijk[1].last() * n_ref, this->ijk[2].last() * n_ref},    // ubound
            {n_ref, n_ref, n_ref}, // stride
            }
        {
          assert(p.n_fra_iter > 0);
          assert(n_ref & 2 == 0);

          for (int d = 0; d < ct_params_t::n_dims; ++d)
          {
            if(d==1)
            {
              ijk_ref.lbound(d) = this->mem->slab(this->mem->grid_size_ref[d], this->rank, this->mem->size).first();// this->ijk.lbound(d) * n_ref;
              ijk_ref.ubound(d) = this->mem->slab(this->mem->grid_size_ref[d], this->rank, this->mem->size).last();
            }
            else
            {
              ijk_ref.lbound(d) = this->mem->grid_size_ref[d].first();
              ijk_ref.ubound(d) = this->mem->grid_size_ref[d].last();
            }
          }
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
