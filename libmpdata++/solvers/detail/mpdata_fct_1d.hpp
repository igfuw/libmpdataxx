/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief FCT option for MPDATA as formulated in @copybrief Smolarkiewicz_and_Grabowski_1990
 */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_fct_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fct_1d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <typename ct_params_t, int minhalo> 
      class mpdata_fct<
	ct_params_t, 
	minhalo,
	typename std::enable_if<ct_params_t::n_dims == 1>::type
      > : public detail::mpdata_fct_common<ct_params_t, minhalo> 
      {
	using parent_t = detail::mpdata_fct_common<ct_params_t, minhalo>; 
	using parent_t::parent_t; // inheriting constructors

	void fct_init(int e)
	{
	  const auto i = this->i^1; // TODO: isn't it a race condition with more than one thread?
	  const auto psi = this->mem->psi[e][this->n[e]];

	  /// \f$ \psi^{max}_{i}=max_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n  
	  /// \f$ \psi^{min}_{i}=min_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n    
	  /// eq.(20a, 20b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
	  this->psi_min(i) = min(min(psi(i-1), psi(i)), psi(i+1));
	  this->psi_max(i) = max(max(psi(i-1), psi(i)), psi(i+1)); 
	}

	void fct_adjust_antidiff(int e, int iter)
	{
	  const int d = 0; // 1D version -> working in x dimension only
	  const auto psi = this->mem->psi[e][this->n[e]];
	  const auto &GC_corr = parent_t::GC_corr(iter);
	  const auto &G = *this->mem->G;
	  const auto &im = this->im; // calculating once for i-1/2 and i+1/2

	  // fill halos in GC_corr
          this->xchng_vctr_alng(GC_corr);

          // calculating betas
          this->beta_up(this->ijk) = formulae::mpdata::beta_up<ct_params_t::opts>(psi, this->psi_max, GC_corr[d], G, this->i);
          this->beta_dn(this->ijk) = formulae::mpdata::beta_dn<ct_params_t::opts>(psi, this->psi_min, GC_corr[d], G, this->i);

          // fill halos for betas
          this->xchng_sclr(this->beta_up);
          this->xchng_sclr(this->beta_dn);

	  // calculating the monotonic corrective velocity
	  this->GC_mono[d]( this->im+h ) = formulae::mpdata::GC_mono<ct_params_t::opts>(
	    psi, this->beta_up, this->beta_dn, GC_corr[d], G, im
	  );
	  
	  // in the last iteration waiting as advop for the next equation will overwrite psi_min/psi_max
	  if (iter == this->n_iters - 1 && parent_t::n_eqns > 1) this->mem->barrier();  // TODO: move to common
	}
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
