/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include <libmpdata++/solvers/detail/boussinesq_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class boussinesq_sgs_common : public boussinesq_common<ct_params_t>
      {
        using parent_t = boussinesq_common<ct_params_t>;

        public:
        using real_t = typename ct_params_t::real_t;

        protected:
        // member fields
        real_t prandtl_num;
        typename parent_t::arr_t &rcdsn_num, &full_tht, &tdef_sq, &mix_len, &hflux_srfc;
        arrvec_t<typename parent_t::arr_t> &grad_tht;

        template <int nd = ct_params_t::n_dims> 
        void calc_rcdsn_num(typename std::enable_if<nd == 2>::type* = 0)
        {
          rcdsn_num(this->ijk) = this->g * 0.5 * (
                                             grad_tht[ct_params_t::n_dims - 1](this->i, this->j - h)
                                           + grad_tht[ct_params_t::n_dims - 1](this->i, this->j + h)
                                           ) / (this->Tht_ref * tdef_sq(this->ijk));
        }

        template <int nd = ct_params_t::n_dims> 
        void calc_rcdsn_num(typename std::enable_if<nd == 3>::type* = 0)
        {
          rcdsn_num(this->ijk) = this->g * 0.5 * (
                                             grad_tht[ct_params_t::n_dims - 1](this->i, this->j, this->k - h)
                                           + grad_tht[ct_params_t::n_dims - 1](this->i, this->j, this->k + h)
                                           ) / (this->Tht_ref * tdef_sq(this->ijk));
        }


        void multiply_sgs_visc()
        {
          static_assert(static_cast<stress_diff_t>(ct_params_t::stress_diff) == compact,
                        "boussinesq with smagorinsky model requires compact sgs formulation");

          // surface indices
          auto ij = this->ijk;
          ij.lbound(ct_params_t::n_dims - 1) = 0;
          ij.ubound(ct_params_t::n_dims - 1) = 0;

          if (!this->mem->G)
          {
            this->hflux_srfc(ij) = -this->hflux_const;
          }
          else
          {
            this->hflux_srfc(ij) = -this->hflux_const * (*this->mem->G)(ij);
          }

          this->calc_full_tht(full_tht);

          this->xchng_pres(full_tht, this->ijk);
          formulae::nabla::calc_grad_cmpct<parent_t::n_dims>(grad_tht, this->full_tht, this->ijk, this->ijkm, this->dijk);
          
          tdef_sq(this->ijk) = formulae::stress::calc_tdef_sq_cmpct<ct_params_t::n_dims>(this->tau, this->ijk);

          calc_rcdsn_num();

          this->k_m(this->ijk) = where(
                                       rcdsn_num(this->ijk) / prandtl_num < 1,
                                       pow(this->smg_c * mix_len(this->ijk), 2)
                                       * sqrt(tdef_sq(this->ijk) * (1 - rcdsn_num(this->ijk) / prandtl_num)),
                                       0
                                      );
          // one level above surface
          auto ijp1 = this->ijk;
          ijp1.lbound(ct_params_t::n_dims - 1) = 1;
          ijp1.ubound(ct_params_t::n_dims - 1) = 1;

          this->k_m(ij) = this->k_m(ijp1);
          
          this->xchng_sclr(this->k_m, this->ijk, 1);

          // havo to use modified ijkm due to shared-memory parallelisation, otherwise overlapping ranges
          // would lead to double multiplications
          // TODO: better way ?
          auto ijkm_aux = this->ijkm;
          if (this->rank > 0)
            ijkm_aux[0] = this->ijk[0];

          formulae::stress::multiply_tnsr_cmpct<ct_params_t::n_dims, ct_params_t::opts>(this->tau, 1.0, this->k_m, *this->mem->G, ijkm_aux);

          this->xchng_sgs_tnsr_offdiag(this->tau, this->tau_srfc, this->ijk, this->ijkm);

          formulae::stress::multiply_vctr_cmpct<ct_params_t::n_dims, ct_params_t::opts>(grad_tht,
                                                                                        1.0 / prandtl_num,
                                                                                        this->k_m,
                                                                                        *this->mem->G,
                                                                                        this->ijk);

          this->xchng_sgs_vctr(grad_tht, hflux_srfc, this->ijk);
          // hack, convinient place to update the heat flux forcing
          this->hflux_frc(this->ijk) = formulae::stress::flux_div_cmpct<parent_t::n_dims, ct_params_t::opts>(grad_tht,
                                                                                                             *this->mem->G,
                                                                                                             this->ijk,
                                                                                                             this->dijk);
        }

        public:
        struct rt_params_t : parent_t::rt_params_t 
        { 
          real_t prandtl_num = 0;
        };

        // ctor
        boussinesq_sgs_common( 
          typename parent_t::ctor_args_t args, 
          const rt_params_t &p
        ) :
          parent_t(args, p),
          prandtl_num(p.prandtl_num),
          rcdsn_num(args.mem->tmp[__FILE__][0][0]),
          full_tht(args.mem->tmp[__FILE__][1][0]),
          tdef_sq(args.mem->tmp[__FILE__][2][0]),
          mix_len(args.mem->tmp[__FILE__][3][0]),
          grad_tht(args.mem->tmp[__FILE__][4]),
          hflux_srfc(args.mem->tmp[__FILE__][5][0])
        {
          //assert(prandtl_num != 0);
        }

        static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
        {
          parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1); // rcdsn_num
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1); // full_tht
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1); // tdef_sq
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1, "mix_len");
          parent_t::alloc_tmp_vctr(mem, __FILE__); // grad_tht
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1, "", true); // hflux_srfc
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
