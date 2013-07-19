/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief FCT option for MPDATA as formulated in @copybrief Smolarkiewicz_and_Grabowski_1990
 */

#pragma once

#include <libmpdata++/solvers/adv/mpdata_2d.hpp>
#include <libmpdata++/solvers/adv/detail/mpdata_fct_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fct_2d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    template <typename real_t, int n_iters, int n_eqs, formulae::mpdata::opts_t opts, int minhalo> 
    class mpdata_fct<real_t, n_iters, 2, n_eqs, opts, minhalo> : 
      public detail::mpdata_fct_common<real_t, n_iters, 2, n_eqs, opts, minhalo> 
    {
      using parent_t = detail::mpdata_fct_common<real_t, n_iters, 2, n_eqs, opts, minhalo>; 

      void fct_init(int e)
      {
	const rng_t i = this->i^1, j = this->j^1; // not optimal - with multiple threads some indices are repeated among threads
	const typename parent_t::arr_t psi = this->state(e); // TODO:! powinno być psi/rho!

	this->psi_min(i,j) = min(min(min(min(
                       psi(i,j+1),
          psi(i-1,j)), psi(i,j  )), psi(i+1,j)),
                       psi(i,j-1)
        );
	this->psi_max(i,j) = max(max(max(max(
                       psi(i,j+1),
          psi(i-1,j)), psi(i,j  )), psi(i+1,j)), 
                       psi(i,j-1)
        ); 
      }

      void fct_adjust_antidiff(int e, int iter)
      {
	const auto &C_corr = parent_t::C_corr(iter);
	const auto psi = this->state(e);
	const auto &im = this->im; // calculating once for i-1/2 and i+1/2
	const auto &jm = this->jm; // calculating once for i-1/2 and i+1/2

	// fill halos -> mpdata works with halo=1, we need halo=2
// TODO: other option would be to define im as a function of halo in mpdata!
	this->mem->barrier();
	this->bcxl->fill_halos_vctr(C_corr[0], this->j^1); // TODO: one xchng call?
	this->bcxr->fill_halos_vctr(C_corr[0], this->j^1);
	this->bcyl->fill_halos_vctr(C_corr[1], this->i^1); // TODO: one xchng call?
	this->bcyr->fill_halos_vctr(C_corr[1], this->i^1);
	this->mem->barrier();

	// calculating the monotonic corrective velocity
	this->C_mono[0]( im+h, jm ) = formulae::mpdata::C_mono<opts, 0>(psi, this->psi_min, this->psi_max, C_corr, im, jm);
	this->C_mono[1]( im, jm+h ) = formulae::mpdata::C_mono<opts, 1>(psi, this->psi_min, this->psi_max, C_corr, jm, im);
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
