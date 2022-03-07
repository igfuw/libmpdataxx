/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

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
      class mpdata_rhs_vip_prs_sgs_fra : public mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
      {
        using parent_t = mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>;

        public:

        using real_t = typename ct_params_t::real_t;

        protected:

//        const int n_fra; // number of fields with fractal reconstruction
//        constexpr int n_rec_cell = pow(2, ct_params::n_fra_rec) ; // number of reconstructed cells in each direction

  //      idx_t<ct_params_t::n_dims> ijkm;

        public:

        struct rt_params_t : parent_t::rt_params_t
        {
          unsigned int n_fra_iter = 1; // number of iterations of fractal reconstruction
        };

        // ctor
        mpdata_rhs_vip_prs_sgs_fra(
          typename parent_t::ctor_args_t args,
          const rt_params_t &p
        ) :
          parent_t(args, p)
        {
        //  for (int d = 0; d < ct_params_t::n_dims; ++d)
        }

        static void alloc(
          typename parent_t::mem_t *mem,
          const int &n_iters
        ) {
          parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, ct_params_t::n_eqns, mem->grid_size_ref, "", false); // rec_psi
          parent_t::alloc_tmp_vctr(mem, __FILE__, mem->grid_size_ref);                      // rec_GC

          mem->psi_ref = mem->tmp[__FILE__][0];
          mem->GC_ref = mem->tmp[__FILE__][1];
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
