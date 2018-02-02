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
	  const auto i1 = this->i^1; // TODO: isn't it a race condition with more than one thread?
	  const auto psi = this->mem->psi[e][this->n[e]];

	  /// \f$ \psi^{max}_{i}=max_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n  
	  /// \f$ \psi^{min}_{i}=min_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n    
	  /// eq.(20a, 20b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
	  this->psi_min(i1) = min(min(psi(i1-1), psi(i1)), psi(i1+1));
	  this->psi_max(i1) = max(max(psi(i1-1), psi(i1)), psi(i1+1)); 
	}

	void fct_adjust_antidiff(int e, int iter)
	{
	  const int d = 0; // 1D version -> working in x dimension only
	  const auto psi = this->mem->psi[e][this->n[e]];
	  auto &GC_corr = parent_t::GC_corr(iter);
	  const auto &G = *this->mem->G;
	  const auto &im = this->im; // calculating once for i-1/2 and i+1/2
	  const auto i1 = this->i^1; 
	  const auto im1 = this->im^1; 

	  // fill halos in GC_corr
          this->xchng_vctr_alng(GC_corr, true);

          // calculation of fluxes for betas denominators
          if (opts::isset(ct_params_t::opts, opts::iga))
          {
            this->flux_ptr = &GC_corr;
          }
          else
          {
            this->flux[0](im1+h) = formulae::donorcell::make_flux<ct_params_t::opts>(psi, GC_corr[0], im1);
            this->flux_ptr = &this->flux;
          }

          const auto &flx = (*(this->flux_ptr));

          // sanity check for input
          assert(std::isfinite(sum(flx[0](i1^h))));

          // calculating betas
          formulae::mpdata::beta_up<ct_params_t::opts>(this->beta_up, psi, this->psi_max, flx, G, i1);
          formulae::mpdata::beta_dn<ct_params_t::opts>(this->beta_dn, psi, this->psi_min, flx, G, i1);
	  
          // assuring flx, psi_min and psi_max are not overwritten
          this->beta_barrier(iter);

	  // calculating the monotonic corrective velocity
	  formulae::mpdata::GC_mono<ct_params_t::opts>(this->GC_mono[d], psi, this->beta_up, this->beta_dn, GC_corr[d], G, im);
	}
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
