/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_fra_3d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    struct mpdata_rhs_vip_prs_sgs_fra_family_tag {};

    template<typename ct_params_t, int minhalo = 0, class enableif = void>
    class mpdata_rhs_vip_prs_sgs_fra;

    // no fields with fractal reconstruction, inherit from mpdata_rhs_vip_prs_sgs
    template<typename ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs_sgs_fra<
      ct_params_t, minhalo,
      typename std::enable_if_t<(int)ct_params_t::fractal_recon == (int)0> 
    > : public mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
    {
      using parent_t = mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors
    };

    template <class ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs_sgs_fra<
      ct_params_t, minhalo,
      typename std::enable_if_t<(int)ct_params_t::fractal_recon != (int)0>
    > : public detail::mpdata_rhs_vip_prs_sgs_fra_dim<ct_params_t, minhalo>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_sgs_fra_dim<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors

      protected:
      using solver_family = mpdata_rhs_vip_prs_sgs_fra_family_tag;

      static void alloc(
        typename parent_t::mem_t *mem,
        const int &n_iters
      ) {
        parent_t::alloc(mem, n_iters);
        parent_t::alloc_tmp_sclr_ref(mem, __FILE__, ct_params_t::n_eqns); // psi_ref
//        parent_t::alloc_tmp_vctr(mem, __FILE__, mem->grid_size_ref);                                 // GC_ref

        mem->psi_ref = mem->tmp[__FILE__][0];
//        mem->GC_ref = mem->tmp[__FILE__][1];
      }
    };
  } // namespace solvers
} // namespace libmpdataxx
