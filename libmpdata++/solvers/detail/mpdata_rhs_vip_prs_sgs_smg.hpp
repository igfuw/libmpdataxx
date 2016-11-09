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
      template <class ct_params_t>
      class mpdata_rhs_vip_prs_sgs_smg : public detail::mpdata_rhs_vip_prs_sgs_common<ct_params_t>
      {
	using parent_t = detail::mpdata_rhs_vip_prs_sgs_common<ct_params_t>;

        protected:
        
	typename parent_t::real_t smg_c;
        typename parent_t::arr_t &k_m;

        void multiply_sgs_visc()
        {
          auto dlta = std::accumulate(this->dijk.begin(), this->dijk.end(), 0.) / 3;
          k_m(this->ijk) = pow(smg_c * dlta, 2) *  formulae::stress::calc_tdef_sq<ct_params_t::n_dims>(this->tau, this->ijk);

          for (auto& t : this->tau)
          {
            t(this->ijk) *= k_m(this->ijk);
          }
        }

	public:

        struct rt_params_t : parent_t::rt_params_t
        {
          typename ct_params_t::real_t smg_c = 0.15;
        };

	// ctor
	mpdata_rhs_vip_prs_sgs_smg(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
          smg_c(p.smg_c),
          k_m(args.mem->tmp[__FILE__][0][0])
	{
          if (smg_c == 0) throw std::runtime_error("smg_c == 0");
        }

        static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
          parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1); // k_m
        }
      }; 
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
