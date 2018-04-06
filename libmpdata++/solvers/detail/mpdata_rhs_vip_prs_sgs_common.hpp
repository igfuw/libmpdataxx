/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#pragma once

#include <numeric>
#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/formulae/idxperm.hpp>
#include <libmpdata++/formulae/stress_formulae.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    enum stress_diff_t
    {
      normal,
      pade,
      compact
    };

    const std::map<stress_diff_t, std::string> sdiff2string = {
      {normal , "normal" },
      {pade   , "pade"   },
      {compact, "compact"}
    };

    namespace detail
    {
      template <class ct_params_t, int minhalo>
      class mpdata_rhs_vip_prs_sgs_common : public mpdata_rhs_vip_prs<ct_params_t, minhalo>
      {
        using parent_t = mpdata_rhs_vip_prs<ct_params_t, minhalo>;

        protected:

        // member fields
        arrvec_t<typename parent_t::arr_t> &tau;
        arrvec_t<typename parent_t::arr_t> &tau_srfc;
        typename parent_t::arr_t           &vip_div;
        arrvec_t<typename parent_t::arr_t> &drv;
        arrvec_t<typename parent_t::arr_t> &wrk;

        std::array<rng_t, ct_params_t::n_dims> ijkm;
        typename parent_t::real_t cdrag;

        virtual void multiply_sgs_visc() = 0;

        // apply pade scheme to the d-th element of the drv array
        void pade_deriv(int d)
        {
          wrk[0](this->ijk) = drv[d](this->ijk);
          for (int m = 0; m < 3; ++m)
          {
            this->xchng_sclr(wrk[0], this->ijk);
            formulae::stress::pade_dispatch<ct_params_t::n_dims>(wrk, this->ijk, d);
            // finish calculation of wrk[1] before modyfying wrk[0]
            this->mem->barrier();
            wrk[0](this->ijk) += (drv[d](this->ijk) - 0.25 * wrk[1](this->ijk));
          }
          drv[d](this->ijk) = wrk[0](this->ijk);
          // needed because otherwise other threads could start calculating pade correction
          // to the next derivative
          this->mem->barrier();
        }

        void vip_rhs_expl_calc()
        {
          parent_t::vip_rhs_expl_calc();

          using ix = typename ct_params_t::ix;
          using namespace arakawa_c;

          // TODO: get rid of superfluous barriers
          for (auto& vip : this->vips())
            this->xchng_sclr(vip, this->ijk, 1);
          
          if (static_cast<stress_diff_t>(ct_params_t::stress_diff) == compact)
          {
            formulae::stress::calc_drag_cmpct<ct_params_t::n_dims, ct_params_t::opts>(tau_srfc,
                                                                                      this->vips(),
                                                                                      *this->mem->G,
                                                                                      cdrag,
                                                                                      this->ijk,
                                                                                      ijkm);

            if (this->mem->G)
            {
              formulae::stress::calc_vip_div_cmpct<ct_params_t::n_dims>(vip_div, this->vips(), *this->mem->G, this->ijk, this->dijk);
              this->xchng_sgs_div(vip_div, this->ijk);
            }
            else
            {
              // TODO: do not use vip_div when G == 1
              vip_div(this->ijk) = 0;
              this->xchng_sgs_div(vip_div, this->ijk);
            }

            formulae::stress::calc_deform_cmpct<ct_params_t::n_dims>(tau, this->vips(), vip_div, this->ijk, ijkm, this->dijk);

            this->xchng_sgs_tnsr_diag(tau, this->vips()[ct_params_t::n_dims - 1], vip_div, this->ijk);
            this->xchng_sgs_tnsr_offdiag(tau, tau_srfc, this->ijk, this->ijkm);
            
            // multiply deformation tensor by sgs viscosity to obtain stress tensor
            multiply_sgs_visc();
            
            // update forces
            formulae::stress::calc_stress_rhs_cmpct<ct_params_t::n_dims, ct_params_t::opts>(this->vip_rhs,
                                                                                            tau,
                                                                                            *this->mem->G,
                                                                                            this->ijk,
                                                                                            this->dijk,
                                                                                            2.0);
          }
          else
          {
            // calculate velocity gradient tensor
            formulae::stress::calc_vgrad<ct_params_t::n_dims>(drv, this->vips(), this->ijk, this->dijk);

            // optionally correct derivatives using Pade scheme
            if ((stress_diff_t)ct_params_t::stress_diff == pade)
            {
              for (int d = 0; d < std::pow(static_cast<int>(ct_params_t::n_dims), 2); ++d)
                pade_deriv(d);
            }
            
            // calculate independent components of deformation tensor
            formulae::stress::calc_deform<ct_params_t::n_dims>(tau, drv, this->ijk);
            
            // multiply deformation tensor by sgs viscosity to obtain stress tensor
            multiply_sgs_visc();
            
            // TODO: get rid of superfluous barriers
            for (auto& t : tau)
            {
              this->xchng_sclr(t, this->ijk);
            }
            // calculate elements of stress tensor divergence
            formulae::stress::calc_stress_div<ct_params_t::n_dims>(drv, tau, this->ijk, this->dijk);

            // optionally correct derivatives using Pade scheme
            if ((stress_diff_t)ct_params_t::stress_diff == pade)
            {
              for (int d = 0; d < std::pow(static_cast<int>(ct_params_t::n_dims), 2); ++d)
                pade_deriv(d);
            }

            // update forces
            formulae::stress::calc_stress_rhs<ct_params_t::n_dims>(this->vip_rhs, drv, this->ijk, 2.0);
          }
        }

        public:

        struct rt_params_t : parent_t::rt_params_t 
        { 
          typename parent_t::real_t cdrag = 0; 
        };

        // ctor
        mpdata_rhs_vip_prs_sgs_common(
          typename parent_t::ctor_args_t args,
          const rt_params_t &p
        ) :
          parent_t(args, p),
          tau(args.mem->tmp[__FILE__][0]),
          tau_srfc(args.mem->tmp[__FILE__][1]),
          vip_div(args.mem->tmp[__FILE__][2][0]),
          drv(args.mem->tmp[__FILE__][3]),
          wrk(args.mem->tmp[__FILE__][4]),
          cdrag(p.cdrag)
        {
          for (int d = 0; d < ct_params_t::n_dims; ++d)
          {
            ijkm[d] = rng_t(this->ijk[d].first() - 1, this->ijk[d].last());
          }
        }

        static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
          parent_t::alloc(mem, n_iters);
          // no staggering for non-compact differencing  
          if (static_cast<stress_diff_t>(ct_params_t::stress_diff) != compact)
          {
            parent_t::alloc_tmp_sclr(mem, __FILE__, 3 * (ct_params_t::n_dims - 1)); // unique strain rate tensor elements
            parent_t::alloc_tmp_sclr(mem, __FILE__, ct_params_t::n_dims, "", true); // unstaggered tau_srfc
          }
          else
          {
            if (ct_params_t::n_dims == 2)
            {
              parent_t::alloc_tmp_stgr(mem,
                                       __FILE__,
                                       3, // unique strain rate tensor elements
                                       {{true, false}, {false, true}, {true, true}}
                                      );
              parent_t::alloc_tmp_stgr(mem, __FILE__, 1, {{true, false}}, true); // tau_srfc
            }
            else
            {
              parent_t::alloc_tmp_stgr(mem,
                                       __FILE__,
                                       6, // unique strain rate tensor elements
                                       {{true, false, false}, {false, true, false}, {false, false, true},
                                        {true, true, false}, {true, false, true}, {false, true, true}}
                                      );
              parent_t::alloc_tmp_stgr(mem, __FILE__, 2, {{true, false, false}, {false, true, false}}, true); // tau_srfc
            }
          }
          if (ct_params_t::n_dims == 2)
          {
            parent_t::alloc_tmp_stgr(mem, __FILE__, 1, {{false, true}}); // vip_div
          }
          else
          {
            parent_t::alloc_tmp_stgr(mem, __FILE__, 1, {{false, false, true}}); // vip_div
          }
          // TODO: do not allocate unnecessary memory when not using pade differencing
          parent_t::alloc_tmp_sclr(mem, __FILE__, std::pow(static_cast<int>(ct_params_t::n_dims), 2)); // drv
          parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // wrk
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
