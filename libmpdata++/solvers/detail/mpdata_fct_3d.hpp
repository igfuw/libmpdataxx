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
	  const auto i1 = this->i^1, j1 = this->j^1, k1 = this->k^1; // not optimal - with multiple threads some indices are repeated among threads
	  const auto psi = this->mem->psi[e][this->n[e]]; 

	  this->psi_min(i1,j1,k1) = min(min(min(min(min(min(
			psi(i1,  j1,  k1),
			psi(i1+1,j1,  k1)),
			psi(i1-1,j1,  k1)),
			psi(i1,  j1+1,k1)),
			psi(i1,  j1-1,k1)),
			psi(i1,  j1,  k1+1)),
			psi(i1,  j1,  k1-1)
	  );
			
	  this->psi_max(i1,j1,k1) = max(max(max(max(max(max(
			psi(i1,  j1,  k1),
			psi(i1+1,j1,  k1)), 
			psi(i1-1,j1,  k1)),
			psi(i1,  j1+1,k1)),
			psi(i1,  j1-1,k1)),
			psi(i1,  j1,  k1+1)), 
			psi(i1,  j1,  k1-1) 
	  ); 
	}

	void fct_adjust_antidiff(int e, int iter)
	{
	  const auto psi = this->mem->psi[e][this->n[e]];
	  auto &GC_corr = parent_t::GC_corr(iter);
          const auto &G = *this->mem->G;
	  const auto &im(this->im), &jm(this->jm), &km(this->km); // calculating once for (i/j/k)-1/2 and (i/j/k)+1/2

          // not optimal - with multiple threads some indices are repeated among threads
	  const auto 
            i = this->i, j = this->j, k = this->k,
            i1 = i^1, j1 = j^1, k1 = k^1,
            im1 = this->im^1, jm1 = this->jm^1, km1 = this->km^1; 

	  // fill halos -> mpdata works with halo=1, we need halo=2
          this->xchng_vctr_alng(GC_corr, true);
          this->xchng_vctr_nrml(this->GC_corr(iter), this->ijk);
          
          // calculation of fluxes for betas denominators
          if (opts::isset(ct_params_t::opts, opts::iga))
          {
            this->flux_ptr = &GC_corr;
          }
          else
          {
            this->flux[0](im1+h, j1,    k1   ) = formulae::donorcell::make_flux<ct_params_t::opts, 0>(psi, GC_corr[0], im1, j1, k1);
            this->flux[1](i1,    jm1+h, k1   ) = formulae::donorcell::make_flux<ct_params_t::opts, 1>(psi, GC_corr[1], jm1, k1, i1);
            this->flux[2](i1,    j1,    km1+h) = formulae::donorcell::make_flux<ct_params_t::opts, 2>(psi, GC_corr[2], km1, i1, j1);
            this->flux_ptr = &this->flux;
          }


          const auto &flx = (*(this->flux_ptr));

          // calculating betas
          formulae::mpdata::beta_up<ct_params_t::opts>(this->beta_up, psi, this->psi_max, flx, G, i1, j1, k1);
          formulae::mpdata::beta_dn<ct_params_t::opts>(this->beta_dn, psi, this->psi_min, flx, G, i1, j1, k1);
        

          // should detect the need for ext=1 in hallo-filling above
	  assert(std::isfinite(sum(this->beta_up(i1, j, k))));
	  assert(std::isfinite(sum(this->beta_up(i, j1, k))));
	  assert(std::isfinite(sum(this->beta_up(i, j, k1))));
	  assert(std::isfinite(sum(this->beta_dn(i1, j, k))));
	  assert(std::isfinite(sum(this->beta_dn(i, j1, k))));
	  assert(std::isfinite(sum(this->beta_dn(i, j, k1))));

          // assuring flx, psi_min and psi_max are not overwritten
          this->beta_barrier(iter);

	  // calculating the monotonic corrective velocity
	  formulae::mpdata::GC_mono<ct_params_t::opts, 0>(this->GC_mono, psi, this->beta_up, this->beta_dn, GC_corr, G, im, this->j, this->k);
	  formulae::mpdata::GC_mono<ct_params_t::opts, 1>(this->GC_mono, psi, this->beta_up, this->beta_dn, GC_corr, G, jm, this->k, this->i);
	  formulae::mpdata::GC_mono<ct_params_t::opts, 2>(this->GC_mono, psi, this->beta_up, this->beta_dn, GC_corr, G, km, this->i, this->j);
        }

      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
