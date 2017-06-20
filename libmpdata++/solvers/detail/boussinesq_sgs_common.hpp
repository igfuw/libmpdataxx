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
        real_t prandtl_num, hflux_surfc;
        typename parent_t::arr_t &rcdsn_num, &full_tht, &tdef_sq, &mix_len;
        arrvec_t<typename parent_t::arr_t> &grad_tht;

        void multiply_sgs_visc()
        {
          this->calc_full_tht(full_tht);

          this->xchng_sclr(full_tht, this->ijk);
          formulae::nabla::calc_grad<parent_t::n_dims>(grad_tht, this->full_tht, this->ijk, this->dijk);
          
          tdef_sq(this->ijk) = formulae::stress::calc_tdef_sq<ct_params_t::n_dims>(this->tau, this->ijk);

          rcdsn_num(this->ijk) = this->g * grad_tht[ct_params_t::n_dims - 1](this->ijk) / (this->Tht_ref * tdef_sq(this->ijk));

          //auto dlta = std::accumulate(this->dijk.begin(), this->dijk.end(), 0.) / 3;

          this->k_m(this->ijk) = where(
                                       rcdsn_num(this->ijk) / prandtl_num < 0,
                                       pow(this->smg_c * mix_len(this->ijk), 2)
                                       * sqrt(tdef_sq(this->ijk) * (1 - rcdsn_num(this->ijk) / prandtl_num)),
                                       0
                                      );

          for (auto& t : this->tau)
          {
            t(this->ijk) *= this->k_m(this->ijk);
          }
         
          for (auto& g : grad_tht)
          {
            g(this->ijk) *= -this->k_m(this->ijk) / prandtl_num;
          }

          for (int d = 0; d < ct_params_t::n_dims; ++d)
          {
            // surface indices
            auto ij = this->ijk;
            ij.lbound(ct_params_t::n_dims - 1) = 0;
            ij.ubound(ct_params_t::n_dims - 1) = 0;
            grad_tht[d](ij) = (d == ct_params_t::n_dims - 1) ? hflux_surfc : 0;
            this->xchng_sclr(grad_tht[d], this->ijk);
          }

          // hack, convinient place to update the heat_flux
          this->hflux(this->ijk) = formulae::nabla::div<parent_t::n_dims>(grad_tht, this->ijk, this->dijk);
        }

        //template<bool is_smg = (ct_params_t::sgs_scheme == smg)>
        //void multiply_sgs_visc_impl(typename std::enable_if<is_smg>::type* = 0)
        //{
        //}
        //
        //template<bool is_smg = (ct_params_t::sgs_scheme == smg)>
        //void multiply_sgs_visc_impl(typename std::enable_if<!is_smg>::type* = 0) {}

        public:
        struct rt_params_t : parent_t::rt_params_t 
        { 
          real_t prandtl_num = 0, hflux_surfc = 0; 
        };

        // ctor
        boussinesq_sgs_common( 
          typename parent_t::ctor_args_t args, 
          const rt_params_t &p
        ) :
          parent_t(args, p),
          prandtl_num(p.prandtl_num),
          hflux_surfc(p.hflux_surfc),
          rcdsn_num(args.mem->tmp[__FILE__][0][0]),
          full_tht(args.mem->tmp[__FILE__][1][0]),
          tdef_sq(args.mem->tmp[__FILE__][2][0]),
          mix_len(args.mem->tmp[__FILE__][3][0]),
          grad_tht(args.mem->tmp[__FILE__][4])
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
          parent_t::alloc_tmp_sclr(mem, __FILE__, ct_params_t::n_dims); // grad_tht
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
