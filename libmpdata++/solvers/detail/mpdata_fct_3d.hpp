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
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fct_3d.hpp>

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
	typename std::enable_if<ct_params_t::n_dims == 3>::type
      > : public detail::mpdata_fct_common<ct_params_t, minhalo> 
      {
	using parent_t = detail::mpdata_fct_common<ct_params_t, minhalo>; 
	using parent_t::parent_t; // inheriting ctors

	void fct_init(int e)
	{
	  const auto i = this->i^1, j = this->j^1, k = this->k^1; // not optimal - with multiple threads some indices are repeated among threads
	  const auto psi = this->mem->psi[e][this->n[e]]; 

	  this->psi_min(i,j,k) = min(min(min(min(min(min(
			psi(i,j,k),
			psi(i+1,j,k)),
			psi(i-1,j,k)),
			psi(i,j+1,k)),
			psi(i,j-1,k)),
			psi(i,j,k+1)),
			psi(i,j,k-1)
	  );
			
	  this->psi_max(i,j,k) = max(max(max(max(max(max(
			psi(i,j,k),
			psi(i+1,j,k)), 
			psi(i-1,j,k)),
			psi(i,j+1,k)),
			psi(i,j-1,k)),
			psi(i,j,k+1)), 
			psi(i,j,k-1) 
	  ); 
	}

	void fct_adjust_antidiff(int e, int iter)
	{
	  const auto &GC_corr = parent_t::GC_corr(iter);
	  const auto psi = this->mem->psi[e][this->n[e]];
	  const auto &im = this->im; // calculating once for i-1/2 and i+1/2
	  const auto &jm = this->jm; // calculating once for j-1/2 and j+1/2
	  const auto &km = this->km; // calculating once for k-1/2 and k+1/2

	  // fill halos -> mpdata works with halo=1, we need halo=2
  // TODO: other option would be to define im as a function of halo in mpdata!
	  this->mem->barrier();
	  this->bcxl->fill_halos_vctr(GC_corr[0], this->j^1, this->k^1); // TODO: one xchng call?
	  this->bcxr->fill_halos_vctr(GC_corr[0], this->j^1, this->k^1);
	  this->bcyl->fill_halos_vctr(GC_corr[1], this->k^1, this->i^1); // TODO: one xchng call?
	  this->bcyr->fill_halos_vctr(GC_corr[1], this->k^1, this->i^1);
	  this->bczl->fill_halos_vctr(GC_corr[2], this->i^1, this->j^1); // TODO: one xchng call?
	  this->bczr->fill_halos_vctr(GC_corr[2], this->i^1, this->j^1);
	  this->mem->barrier();

	  // calculating the monotonic corrective velocity
	  this->GC_mono[0]( im+h, jm, km ) = formulae::mpdata::C_mono<ct_params_t::opts, 0>(psi, this->psi_min, this->psi_max, GC_corr, im, jm, km);
	  this->GC_mono[1]( im, jm+h, km ) = formulae::mpdata::C_mono<ct_params_t::opts, 1>(psi, this->psi_min, this->psi_max, GC_corr, jm, km, im);
	  this->GC_mono[2]( im, jm, km+h ) = formulae::mpdata::C_mono<ct_params_t::opts, 2>(psi, this->psi_min, this->psi_max, GC_corr, km, im, jm);

	  // in the last iteration waiting as advop for the next equation will overwrite psi_min/psi_max
	  if (iter == this->n_iters - 1) this->mem->barrier();  // TODO: move to common // TODO: different condition in 1D and 2D!
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
