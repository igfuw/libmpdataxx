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

      protected:
      using solver_family = mpdata_rhs_vip_prs_sgs_fra_family_tag;

      //ctor
      mpdata_rhs_vip_prs_sgs_fra(
        typename parent_t::ctor_args_t args,
        const typename parent_t::rt_params_t &p
      ) :
        parent_t(args, p)
      {
        if(args.mem->n_ref > 1) // n_ref==1 doesnt do refinement
        {
          this->c_j.reference(args.mem->tmp[__FILE__][1][0]);
          this->d_j.reference(args.mem->tmp[__FILE__][1][1]);
          this->f_j.reference(args.mem->tmp[__FILE__][1][2]);
        }
      }

      // why alloc here?
      static void alloc(
        typename parent_t::mem_t *mem,
        const int &n_iters
      ) {
        parent_t::alloc(mem, n_iters);

        if(mem->n_ref > 1)
        {
          parent_t::alloc_tmp_sclr_ref(mem, __FILE__, opts::nset(ct_params_t::fractal_recon));
          parent_t::alloc_tmp_sclr_ref(mem, __FILE__, 3); // c_j, d_j, f_j
          mem->psi_ref = mem->tmp[__FILE__][0];
//        parent_t::alloc_tmp_vctr(mem, __FILE__, mem->grid_size_ref);                                 // GC_ref
//        mem->GC_ref = mem->tmp[__FILE__][1];
        }
        else if(mem->n_ref == 1) // no refinement, dont allocate anything, psi_ref is a reference for psi
        {
          for(opts::opts_t e=0; e<ct_params_t::n_eqns; ++e) // loop over all variables
          {
            if(opts::isset(ct_params_t::fractal_recon, opts::bit(e))) // if the variable can be refined
            {
              mem->psi_ref.at(ix_r2r.at(e)).reference(mem->psi.at(e));
            }
          }
        }
      }
    };
  } // namespace solvers
} // namespace libmpdataxx
