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
          formulae::nabla::calc_grad<parent_t::n_dims, true>(grad_tht, this->full_tht, this->ijk, this->dijk);
          
          tdef_sq(this->ijk) = formulae::stress::calc_tdef_sq_cmpct<ct_params_t::n_dims>(this->tau, this->ijk);

          rcdsn_num(this->ijk) = this->g * 0.5 * (
                                             grad_tht[ct_params_t::n_dims - 1](this->i, this->j, this->k - h)
                                           + grad_tht[ct_params_t::n_dims - 1](this->i, this->j, this->k + h)
                                           ) / (this->Tht_ref * tdef_sq(this->ijk));

          //auto dlta = std::accumulate(this->dijk.begin(), this->dijk.end(), 0.) / 3;

          this->k_m(this->ijk) = where(
                                       rcdsn_num(this->ijk) / prandtl_num < 0,
                                       pow(this->smg_c * mix_len(this->ijk), 2)
                                       * sqrt(tdef_sq(this->ijk) * (1 - rcdsn_num(this->ijk) / prandtl_num)),
                                       0
                                      );
          
          this->xchng_sclr(this->k_m, this->ijk);

          const auto &i(this->i), &j(this->j), &k(this->k);
          const auto &im(this->im), &jm(this->jm), &km(this->km);
          const auto &ii = this->rank == 0 ? im : i;

          this->tau[0](ii + h, j, k) *= 0.5 * 
                             (this->k_m(ii + 1, j, k) + this->k_m(ii, j, k));
          
          this->tau[1](ii + h, jm + h, k) *= 0.25 * 
                             (this->k_m(ii + 1, jm, k) + this->k_m(ii, jm, k) +
                             this->k_m(ii + 1, jm + 1, k) + this->k_m(ii, jm + 1, k));
          
          this->tau[2](ii + h, j, km + h) *= 0.25 * 
                             (this->k_m(ii + 1, j, km) + this->k_m(ii, j, km) +
                             this->k_m(ii + 1, j, km + 1) + this->k_m(ii, j, km + 1));
          
          this->tau[3](i, jm + h, k) *= 0.5 * 
                             (this->k_m(i, jm + 1, k) + this->k_m(i, jm, k));
          
          this->tau[4](i, jm + h, km + h) *= 0.25 * 
                             (this->k_m(i, jm, km + 1) + this->k_m(i, jm, km) +
                             this->k_m(i, jm + 1, km + 1) + this->k_m(i, jm + 1, km));

          this->tau[5](i, j, km + h) *= 0.5 * 
                             (this->k_m(i, j, km + 1) + this->k_m(i, j, km));
          //for (auto& t : this->tau)
          //{
          //  t(this->ijk) *= this->k_m(this->ijk);
          //}
         
          //for (auto& g : grad_tht)
          //{
          //  g(this->ijk) *= this->k_m(this->ijk) / prandtl_num;
          //}
          
          grad_tht[0](this->i + h, this->j, this->k) *= 0.5 *
                             (this->k_m(this->i + 1, this->j, this->k) + this->k_m(this->i, this->j, this->k)) /
                             prandtl_num;
          
          grad_tht[1](this->i, this->j + h, this->k) *= 0.5 *
                             (this->k_m(this->i, this->j + 1, this->k) + this->k_m(this->i, this->j, this->k)) /
                             prandtl_num;
          
          grad_tht[2](this->i, this->j, this->k + h) *= 0.5 *
                             (this->k_m(this->i, this->j, this->k + 1) + this->k_m(this->i, this->j, this->k)) /
                             prandtl_num;

          //for (int d = 0; d < ct_params_t::n_dims; ++d)
          //{
          //  // surface indices
          //  auto ij = this->ijk;
          //  ij.lbound(ct_params_t::n_dims - 1) = 0;
          //  ij.ubound(ct_params_t::n_dims - 1) = 0;
          //  grad_tht[d](ij) = (d == ct_params_t::n_dims - 1) ? -hflux_surfc : 0;
          //  this->xchng_sclr(grad_tht[d], this->ijk);
          //}

          this->xchng_vctr_alng(grad_tht, -hflux_surfc);
          // hack, convinient place to update the heat_flux
          this->hflux(this->ijk) = formulae::nabla::div_cmpct<parent_t::n_dims>(grad_tht, this->ijk, this->dijk);
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
          parent_t::alloc_tmp_vctr(mem, __FILE__); // grad_tht
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
