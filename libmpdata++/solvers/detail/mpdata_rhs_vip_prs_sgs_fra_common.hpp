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
        const int n_ref; // number of refinements; refined resolution is dx / n_ref
        
        // refined arrays have no halo, because halo size needs to be known at compile time (why is it so? probably not really necessary?)
//        const int halo_ref; // size of halos of refined scalar (and vector?) arrays

        // TODO: make these const!
        idx_t<ct_params_t::n_dims>  ijk_ref; // range of refinee handled by given solver
        const idxs_t<ct_params_t::n_dims> ijk_r2r; // resolved to refined; refined scalars at the same position as resolved scalars

//        static rng_t rng_sclr_ref(const rng_t &rng) { return rng^halo_ref; }

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
