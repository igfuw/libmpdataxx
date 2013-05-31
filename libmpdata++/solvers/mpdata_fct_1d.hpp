/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief FCT option for MPDATA as formulated in @copybrief Smolarkiewicz_and_Grabowski_1990
 */

#pragma once

#include <libmpdata++/solvers/mpdata_1d.hpp>
#include <libmpdata++/solvers/detail/mpdata_fct_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fct_1d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    template <
      typename real_t, 
      int n_iters, 
      int n_eqs, 
      int halo
    > 
    class mpdata_fct<real_t, n_iters, 1, n_eqs, halo> : public detail::mpdata_fct_common<real_t, n_iters, 1, n_eqs, halo> 
    {
      using parent_t = detail::mpdata_fct_common<real_t, n_iters, 1, n_eqs, halo>; 

      void fct_init(int e)
      {
	const rng_t i = this->i^1; // TODO: isn't it a race condition with more than one thread?
	const typename parent_t::arr_t psi = this->state(e); // TODO:! powinno być psi/rho!

	this->psi_min(i) = min(min(psi(i-1), psi(i)), psi(i+1));
	this->psi_max(i) = max(max(psi(i-1), psi(i)), psi(i+1)); 
      }

      void fct_adjust_antidiff(int e, int iter)
      {
        const int d = 0; // 1D version -> working in x dimension only
        const auto &C_corr = parent_t::C_corr(iter);
        const auto psi = this->state(e); 
        const auto &im = this->im; // calculating once for i-1/2 and i+1/2

        // fill halos -> mpdata works with halo=1, we need halo=2
// TODO: other option would be to define im as a function of halo in mpdata!
        this->mem->barrier();
	this->bcxl->fill_halos_vctr(C_corr[d]); // TODO: one xchng call?
	this->bcxr->fill_halos_vctr(C_corr[d]);
	this->mem->barrier();

        // calculating the monotonic corrective velocity
        this->C_mono[d]( im+h ) = C_corr[d]( im+h ) * where(
          C_corr[d]( im+h ) > 0,
          min(1, min(
            formulae::mpdata::beta_dn(psi, this->psi_min, C_corr[d], im),
            formulae::mpdata::beta_up(psi, this->psi_max, C_corr[d], im + 1)
          )),
          min(1, min(
            formulae::mpdata::beta_up(psi, this->psi_max, C_corr[d], im),
            formulae::mpdata::beta_dn(psi, this->psi_min, C_corr[d], im + 1)
          ))
        );
        // TODO: positive_definite option (z parent_t)
      }

      public:

      // ctor (TODO: C++11 ctor inheritance?)
      mpdata_fct(
        typename parent_t::ctor_args_t args, 
        const typename parent_t::params_t &p
      ) : 
        parent_t(args, p)
      {}   
    };
  }; // namespace solvers
}; // namespace libmpdataxx
