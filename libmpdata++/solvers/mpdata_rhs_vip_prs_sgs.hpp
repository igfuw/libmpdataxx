/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/formulae/idxperm.hpp>
#include <libmpdata++/formulae/stress_formulae.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    enum sgs_scheme_t
    {
      dns,
      smg,
      tke
    };

    enum stress_diff_t
    {
      normal,
      pade
    };

    template <class ct_params_t>
    class mpdata_rhs_vip_prs_sgs : public mpdata_rhs_vip_prs<ct_params_t>
    {
      using parent_t = mpdata_rhs_vip_prs<ct_params_t>;

      protected:

      // member fields
      arrvec_t<typename parent_t::arr_t> &tau;
      arrvec_t<typename parent_t::arr_t> &drv;
      arrvec_t<typename parent_t::arr_t> &wrk;
      typename ct_params_t::real_t eta = 0;

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
        // to the next derivataive
        this->mem->barrier();
      }

      void vip_rhs_expl_calc()
      {
        parent_t::vip_rhs_expl_calc();
        using ix = typename ct_params_t::ix;

        // TODO: get rid of superfluous barriers
        for (auto& vip : this->vips)
          this->xchng_sclr(vip, this->ijk);

        // calculate velocity gradient tensor
        formulae::stress::calc_vgrad<ct_params_t::n_dims>(drv, this->vips, this->ijk, this->dijk);

        // optionally correct derivatives using Pade scheme
        if ((stress_diff_t)ct_params_t::stress_diff == pade)
        {
          for (int d = 0; d < std::pow(static_cast<int>(ct_params_t::n_dims), 2); ++d)
            pade_deriv(d);
        }
        
        // calculate independent components of deformation tensor
        formulae::stress::calc_deform<ct_params_t::n_dims>(tau, drv, this->ijk);

        // TODO: get rid of superfluous barriers
        for (auto& t : tau)
          this->xchng_sclr(t, this->ijk);

        // calculate elements of stress tensor divergence
        formulae::stress::calc_stress_div<ct_params_t::n_dims>(drv, tau, this->ijk, this->dijk);

        // optionally correct derivatives using Pade scheme
        if ((stress_diff_t)ct_params_t::stress_diff == pade)
        {
          for (int d = 0; d < std::pow(static_cast<int>(ct_params_t::n_dims), 2); ++d)
            pade_deriv(d);
        }

        // update forces
        switch ((sgs_scheme_t)ct_params_t::sgs_scheme)
        {
          case dns:
          {
            formulae::stress::calc_stress_rhs<ct_params_t::n_dims>(this->vip_rhs, drv, this->ijk, 2 * eta * this->dt);
            break;
          }
          case smg:
            assert(false);
            break;
          default:
            assert(false);
        }
      }

      public:

      struct rt_params_t : parent_t::rt_params_t
      {
        typename ct_params_t::real_t eta = 0;
      };

      // ctor
      mpdata_rhs_vip_prs_sgs(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) :
	parent_t(args, p),
        eta(p.eta),
        tau(args.mem->tmp[__FILE__][0]),
        drv(args.mem->tmp[__FILE__][1]),
        wrk(args.mem->tmp[__FILE__][2])
      {
        if (eta == 0)
          throw std::runtime_error("eta == 0");
      }

      static void alloc(
        typename parent_t::mem_t *mem, 
        const int &n_iters
      ) {
	parent_t::alloc(mem, n_iters);
        parent_t::alloc_tmp_sclr(mem, __FILE__, 3 * (ct_params_t::n_dims - 1)); // unique strain rate tensor elements
        // TODO: do not allocate unnecessary memory when not using pade differencing
        parent_t::alloc_tmp_sclr(mem, __FILE__, std::pow(static_cast<int>(ct_params_t::n_dims), 2));
        parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // wrk
      }
    };
  } // namespace solvers
} // namespace libmpdataxx
