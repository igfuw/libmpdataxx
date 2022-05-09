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
      // TODO: minimum halo size is 2
      template <class ct_params_t, int minhalo>
      class mpdata_rhs_vip_prs_sgs_fra_common : public mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
      {
        using parent_t = mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>;

        public:

        using real_t = typename ct_params_t::real_t;

        protected:

        typename parent_t::arr_t c_j, d_j, f_j; // parameters used in fractal reconstruction, Akinlabi et al. 2019


        const int n_ref,         // number of refinements; refined resolution is dx / n_ref
                  n_fra_iter;    // number of iterations of grid refinement
        
//        const int halo_ref; // size of halos of refined scalar (and vector?) arrays

        const std::array<int, ct_params_t::n_eqns> ix_r2r;
        constexpr std::array<int, ct_params_t::n_eqns> get_ix_r2r()
        {
          std::array<int, ct_params_t::n_eqns> ix_r2r{};
          int j = 0;
          for(opts::opts_t e=0; e<ct_params_t::n_eqns; ++e)
          {
            ix_r2r.at(e) = opts::isset(ct_params_t::fractal_recon, opts::bit(e)) ? j++ : -1; // -1 index should give segfaults
          }
          return ix_r2r;
        }

        // herlper ranges
        // TODO: make these const!
        idx_t<ct_params_t::n_dims>  ijk_ref, // range of refinee handled by given solver
                                    ijk_ref_with_halo; // same but with a halo in x direction between MPI processes
        const idxs_t<ct_params_t::n_dims> ijk_r2r; // resolved to refined; refined scalars at the same position as resolved scalars

        // range modifying methods used in grid refinement
        // TODO: unify

        static rng_t rng_midpoints(const rng_t &rng, const int rank = 0, const int size = 1) 
        {
          assert(rng.stride() % 2 == 0);
          if(rng.last() - rng.first() < rng.stride()) return rng;
          return rng_t(
            rank == 0 ?        rng.first() + rng.stride() / 2 : rng.first() - rng.stride() / 2,
            rank == size - 1 ? rng.last()  - rng.stride() / 2 : rng.last()  + rng.stride() / 2, 
            rng.stride()); 
        }

        static rng_t rng_midpoints_in_out(const rng_t &rng, const int rank = 0, const int size = 1) 
        {
          assert(rng.stride() % 2 == 0);
          if(rng.last() - rng.first() < rng.stride()) return rng;
          return rng_t(
            rng.first() + rng.stride() / 2,
            rank == size - 1 ? rng.last()  - rng.stride() / 2 : rng.last()  + rng.stride() / 2, 
            rng.stride()); 
        }

        static rng_t rng_midpoints_in(const rng_t &rng, const int rank = 0, const int size = 1) 
        {
          assert(rng.stride() % 4 == 0);
          if(rng.last() - rng.first() < rng.stride()) return rng;
          else return rng_t(
            rank == 0        ? rng.first() - rng.stride() / 4 : rng.first() + rng.stride() / 4,
            rank == size - 1 ? rng.last()  + rng.stride() / 4 : rng.last()  - rng.stride() / 4,
            rng.stride() / 2); 
        }

        static rng_t rng_midpoints_out(const rng_t &rng) 
        {
          assert(rng.stride() % 4 == 0);
          if(rng.last() - rng.first() < rng.stride()) return rng;
          else return rng_t(
            rng.first() - rng.stride() / 4,
            rng.last() + rng.stride() / 4, 
            rng.stride() / 2); 
        }

        static rng_t rng_half_stride(const rng_t &rng, const int rank = 0, const int size = 1) 
        {
          assert(rng.stride() % 2 == 0);
          return rng_t(
            rng.first(), 
            rank < size - 1 ? rng.last() + rng.stride() / 2 : rng.last(), 
            rng.stride() / 2); 
        }

        // reconstruction based on 3 points
        // input is an array pointing to midpoints (to be reconstructed), but with a stride jumping on each one
        // returned range points to the first from the pair of reconstructed points
        // to avoid boundary conditions, in y and z directions it is assumed that the number of points is 3+i*2, i=0,1,2,... (this check is somewhere else)
        // however in the x direction there can be any number of points, because the domain is divided between mpi processes...
        static rng_t rng_dbl_stride(const rng_t &rng) 
        {
          assert(rng.last() != rng.first()); // we need at least 2 midpoints
          assert(rng.stride() % 2 == 0);

          // if domain starts with a point in the middle of a triple - we need to start calculating from a point to the left
          const int offset = ( (rng.first() - rng.stride() / 2) / rng.stride() ) % 2 == 0 ? 0 : - rng.stride();
          const auto first = rng.first() + offset;

          if( ((rng.last() - first) / rng.stride() + 1) % 2 == 0) // even number of midpoints; y and z directions (and sometimes x)
            return rng_t(first, rng.last() - rng.stride(), 2*rng.stride());
          else // odd number of midpoints
            return rng_t(first, rng.last(),                2*rng.stride()); // rely on the halo along x direction
        }

        static rng_t rng_merge(const rng_t &rng1, const rng_t &rng2) 
        {
          assert(rng1.stride() == rng2.stride());
          return rng_t(
            std::min(rng1.first(), rng2.first()),
            std::max(rng1.last(), rng2.last()),
            rng1.stride() / 2);
        }

        // reconstruction based on 3 points, we need up to 2 resolved points between MPI domains
        static rng_t rng_ref_distmem_halo(const rng_t &rng, const int &n_ref, const int rank = 0, const int size = 1)
        {
          return rng_t(
            rank == 0       ? rng.first()                              : rng.first() - (n_ref + n_ref/2),
            rank < size - 1 ? rng.last() + (n_ref + n_ref/2) : rng.last(),
            rng.stride());
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
          n_ref(args.mem->n_ref),
//          halo_ref(n_ref / 2),
          ix_r2r(get_ix_r2r()),
          // ijk_ref init below assumes 3D (shmem decomp dim along y);
          // TODO: move this to concurr_common::init()? add something to ctor_args_t?
          ijk_r2r{
            {this->ijk[0].first() * n_ref, this->ijk[1].first() * n_ref, this->ijk[2].first() * n_ref}, // lbound
            {this->ijk[0].last() * n_ref, this->ijk[1].last() * n_ref, this->ijk[2].last() * n_ref},    // ubound
            {n_ref, n_ref, n_ref}, // stride
            }
        {
          assert(p.n_fra_iter > 0);
          assert(n_ref % 2 == 0);

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
          ijk_ref_with_halo = ijk_ref;
          const auto ijk_ref_with_halo_rng_0 = rng_ref_distmem_halo(ijk_ref[0], this->mem->n_ref, this->mem->distmem.rank(), this->mem->distmem.size());
          ijk_ref_with_halo.lbound(0) = ijk_ref_with_halo_rng_0.first();
          ijk_ref_with_halo.ubound(0) = ijk_ref_with_halo_rng_0.last();
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
