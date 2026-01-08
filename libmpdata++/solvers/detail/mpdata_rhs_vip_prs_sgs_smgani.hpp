/**
  * @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t, int minhalo>
      class mpdata_rhs_vip_prs_sgs_smgani : public detail::mpdata_rhs_vip_prs_sgs_common<ct_params_t, minhalo>
      {
        using parent_t = detail::mpdata_rhs_vip_prs_sgs_common<ct_params_t, minhalo>;

        public:

        using real_t = typename ct_params_t::real_t;

        protected:

        real_t smg_c, c_m;
        arrvec_t<typename parent_t::arr_t> &k_m;

        void multiply_sgs_visc()
        {
          static_assert(ct_params_t::n_dims > 1, "libmpdata++: anisotropic smagorinsky doesn't work in 1D");
          static_assert(static_cast<stress_diff_t>(ct_params_t::stress_diff) == compact, "libmpdata++: anisotropic smagorinsky requires compact stress differencing");

          const auto dlta_h = std::accumulate(this->dijk.begin(), this->dijk.end()-1, real_t(0.)) / (ct_params_t::n_dims-1);
          const auto dlta_v = this->dijk[ct_params_t::n_dims-1];

          // Simon and Chow 2021, eqs. 9 and 10
          k_m[0](this->ijk) = pow(smg_c * dlta_h, 2) * formulae::stress::calc_tdef_sq_cmpct<ct_params_t::n_dims>(this->tau, this->ijk)(this->ijk); // tdef_sq could be cached, but creating a tdef_sq array without similar one in isotropic smg messes with derived classes that do cache tdef_sq also in isotropic... (e.g. boussinesq)
          k_m[1](this->ijk) = pow(smg_c * dlta_v, 2) * formulae::stress::calc_tdef_sq_cmpct<ct_params_t::n_dims>(this->tau, this->ijk)(this->ijk);

          formulae::stress::multiply_tnsr_cmpct<ct_params_t::n_dims, ct_params_t::opts>(this->tau,
                                                                                        real_t(1.0),
                                                                                        this->k_m,
                                                                                        *this->mem->G,
                                                                                        this->ijkm_sep);

          this->xchng_sgs_tnsr_diag(this->tau, this->vips()[ct_params_t::n_dims - 1], this->vip_div, this->ijk);
          this->xchng_sgs_tnsr_offdiag(this->tau, this->tau_srfc, this->ijk, this->ijkm);
        }

        public:

        struct rt_params_t : parent_t::rt_params_t
        {
          real_t smg_c, c_m;
        };

        // ctor
        mpdata_rhs_vip_prs_sgs_smgani(
          typename parent_t::ctor_args_t args,
          const rt_params_t &p
        ) :
          parent_t(args, p),
          smg_c(p.smg_c),
          c_m(p.c_m),
          k_m(args.mem->tmp[__FILE__][0])
        {
          if (smg_c == 0) throw std::runtime_error("libmpdata++: smg_c == 0");
        }

        static void alloc(
          typename parent_t::mem_t *mem,
          const int &n_iters
        ) {
          parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // k_m
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
