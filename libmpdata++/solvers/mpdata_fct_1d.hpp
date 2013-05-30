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
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fct_1d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail 
    {
      // TODO: document why 2
      const int fct_min_halo = 2; // TODO move to fct::formulae?
    }

//    enum rho_enum {rho_constant, rho_profile, rho_variable};

// TODO: mpdata_fct_common
    template <
      typename real_t, 
      int n_iters, 
      int n_eqs = 1, 
//      rho_enum rho_opt = rho_constant, 
      int halo = detail::fct_min_halo
    > 
    class mpdata_fct_1d : public mpdata_1d<real_t, n_iters, n_eqs, detail::max(halo, detail::fct_min_halo)> 
// TODO: inherit (multiply!) from mpdata_fct_common
    {
      public:

      using parent_t = mpdata_1d<real_t, n_iters, n_eqs, detail::max(halo, detail::fct_min_halo)>; 

      static_assert(parent_t::n_iters > 1, "FCT is defined for MPDATA with a corrective iteration (not for donorcell)");

      struct params_t : parent_t::params_t
      {
        // TODO: rho!
      };  

      private:

      typename parent_t::arr_t psi_min, psi_max; // TODO: movo to modata_fct_common
      arrvec_t<typename parent_t::arr_t> C_mono; // TODO: movo to modata_fct_common

      void fct_init(int e)
      {
	const rng_t i = this->i^1; // TODO: isn't it a race condition with more than one thread?
	const typename parent_t::arr_t psi = this->state(e); // TODO:! powinno być psi/rho!

	psi_min(i) = min(min(psi(i-1), psi(i)), psi(i+1));
	psi_max(i) = max(max(psi(i-1), psi(i)), psi(i+1)); 
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
        C_mono[d]( im+h ) = C_corr[d]( im+h ) * where(
          C_corr[d]( im+h ) > 0,
          min(1, min(
            formulae::mpdata::beta_dn(psi, psi_min, C_corr[d], im),
            formulae::mpdata::beta_up(psi, psi_max, C_corr[d], im + 1)
          )),
          min(1, min(
            formulae::mpdata::beta_up(psi, psi_max, C_corr[d], im),
            formulae::mpdata::beta_dn(psi, psi_min, C_corr[d], im + 1)
          ))
        );
        // TODO: positive_definite option (z parent_t)
      }

      arrvec_t<typename parent_t::arr_t> &C(int iter) // TODO: move to mpdata_fct_common
      {
        if (iter > 0) return C_mono;
        return parent_t::C(iter);
      }

      public:

      // ctor
      mpdata_fct_1d(
        typename parent_t::ctor_args_t args, 
        const params_t &p
      ) : 
        parent_t(args, p),
        psi_min(args.mem->tmp[__FILE__][0][0]),
        psi_max(args.mem->tmp[__FILE__][0][1]),
         C_mono(args.mem->tmp[__FILE__][1])
      {}   

      // 1D version
      static void alloc(typename parent_t::mem_t *mem, const int nx)
      {
	parent_t::alloc(mem, nx); 
        parent_t::alloc_tmp_sclr(mem, nx, __FILE__, 2); // psi_min and psi_max
        parent_t::alloc_tmp_vctr(mem, nx, __FILE__);    // C_mono
      }
    };
  }; // namespace solvers
}; // namespace libmpdataxx
