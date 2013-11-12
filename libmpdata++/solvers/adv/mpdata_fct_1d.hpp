/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief FCT option for MPDATA as formulated in @copybrief Smolarkiewicz_and_Grabowski_1990
 */

#pragma once

#include <libmpdata++/solvers/adv/mpdata_1d.hpp>
#include <libmpdata++/solvers/adv/detail/mpdata_fct_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fct_1d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    template <
      typename real_t, 
      int n_iters, 
      int n_eqs, 
      formulae::mpdata::opts_t opts,
      int minhalo
    > 
    class mpdata_fct<real_t, n_iters, 1, n_eqs, opts, minhalo> : 
      public detail::mpdata_fct_common<real_t, n_iters, 1, n_eqs, opts, minhalo> 
    {
      using parent_t = detail::mpdata_fct_common<real_t, n_iters, 1, n_eqs, opts, minhalo>; 

      void fct_init(int e)
      {
	const auto i = this->i^1; // TODO: isn't it a race condition with more than one thread?
	const auto psi = this->state(e);

        /// \f$ \psi^{max}_{i}=max_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n  
        /// \f$ \psi^{min}_{i}=min_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n    
        /// eq.(20a, 20b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
	this->psi_min(i) = min(min(psi(i-1), psi(i)), psi(i+1));
	this->psi_max(i) = max(max(psi(i-1), psi(i)), psi(i+1)); 
      }

      void fct_adjust_antidiff(int e, int iter)
      {
        const int d = 0; // 1D version -> working in x dimension only
        const auto &GC_corr = parent_t::GC_corr(iter);
        const auto psi = this->state(e); 
        const auto &im = this->im; // calculating once for i-1/2 and i+1/2

        // fill halos -> mpdata works with halo=1, we need halo=2
// TODO: other option would be to define im as a function of halo in mpdata!
        this->mem->barrier();
	this->bcxl->fill_halos_vctr(GC_corr[d]); // TODO: one xchng call?
	this->bcxr->fill_halos_vctr(GC_corr[d]);
	this->mem->barrier();

        // calculating the monotonic corrective velocity
	this->GC_mono[d]( im+h ) = formulae::mpdata::C_mono<opts>(psi, this->psi_min, this->psi_max, GC_corr[d], im);
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
