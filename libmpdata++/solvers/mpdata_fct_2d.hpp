/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief FCT option for MPDATA as formulated in @copybrief Smolarkiewicz_and_Grabowski_1990
 */

#pragma once

#include <libmpdata++/solvers/mpdata_2d.hpp>
#include <libmpdata++/solvers/detail/mpdata_fct_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fct_2d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    template <typename real_t, int n_iters, int n_eqs, formulae::mpdata::opts_t opts, int halo> 
    class mpdata_fct<real_t, n_iters, 2, n_eqs, opts, halo> : 
      public detail::mpdata_fct_common<real_t, n_iters, 2, n_eqs, opts, halo> 
    {
      using parent_t = detail::mpdata_fct_common<real_t, n_iters, 2, n_eqs, opts, halo>; 

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
        for (int d = 0; d < parent_t::n_dims; ++d)
	  this->C_mono[d](rng_t::all(), rng_t::all()) = 
	    this->C_corr(iter)[d](rng_t::all(), rng_t::all());
        // TODO
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
