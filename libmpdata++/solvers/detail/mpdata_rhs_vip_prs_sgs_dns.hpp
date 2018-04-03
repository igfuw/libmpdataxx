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
      class mpdata_rhs_vip_prs_sgs_dns : public detail::mpdata_rhs_vip_prs_sgs_common<ct_params_t, minhalo>
      {
	using parent_t = detail::mpdata_rhs_vip_prs_sgs_common<ct_params_t, minhalo>;

        protected:
        
	typename parent_t::real_t eta;

        void multiply_sgs_visc() final
        {
          if (static_cast<stress_diff_t>(ct_params_t::stress_diff) == compact)
          {
            formulae::stress::multiply_tnsr_cmpct<ct_params_t::n_dims>(this->tau, eta, this->ijk);

            this->xchng_sgs_tnsr_diag(this->tau, this->vips()[ct_params_t::n_dims - 1], this->vip_div, this->ijk);
            this->xchng_sgs_tnsr_offdiag(this->tau, this->tau_srfc, this->ijk, this->ijkm);
          }
          else
          {
            for (auto& t : this->tau)
            {
              t(this->ijk) *= eta;
            }
          }
        }

	public:

        struct rt_params_t : parent_t::rt_params_t
        {
          typename ct_params_t::real_t eta = 0;
        };

	// ctor
	mpdata_rhs_vip_prs_sgs_dns(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
          eta(p.eta)
	{
          if (eta == 0) throw std::runtime_error("eta == 0");
        }
      }; 
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
