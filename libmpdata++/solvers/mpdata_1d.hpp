/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_1d.hpp>
#include <libmpdata++/formulae/donorcell_formulae.hpp>
#include <libmpdata++/solvers/detail/solver_1d.hpp>

#include <array>

// TODO: an mpdata_common class?

namespace libmpdataxx
{
  namespace solvers
  {
    using namespace libmpdataxx::arakawa_c;

    template<typename real_t, int n_iters, int n_eqs = 1, int halo = formulae::mpdata::halo>
    class mpdata_1d : public detail::solver_1d<
      real_t, 
      1,
      n_eqs,
      formulae::mpdata::n_tlev,
      detail::max(halo, formulae::mpdata::halo)
    >
    {
      static_assert(n_iters > 0, "n_iters <= 0");

      private: 
      using parent_t = detail::solver_1d< 
        real_t,
        1,
        n_eqs,
        formulae::mpdata::n_tlev,
        detail::max(halo, formulae::mpdata::halo)
      >;

      using C_t = arrvec_t<typename parent_t::arr_t>;

      static const int n_tmp = n_iters > 2 ? 2 : 1; // TODO: this should be in mpdata_common

      // member fields
      std::array<C_t*, n_tmp> tmp;
      rng_t im;

      protected:

      // for Flux-Corrected Transport (TODO: more general names?)
      virtual void fct_init(int e) { }
      virtual void fct_adjust_antidiff(int e) { }

      C_t &C_unco(int iter)
      {
        return (iter == 1) 
	  ? this->mem->C 
	  : (iter % 2) 
	    ? *tmp[1]  // odd iters
	    : *tmp[0]; // even iters
      }

      C_t &C_corr(int iter) 
      {
        return (iter  % 2) 
	  ? *tmp[0]    // odd iters
	  : *tmp[1];   // even iters
      }

      virtual C_t &C(int iter)
      {
        if (iter == 0) return this->mem->C;
        return C_corr(iter);
      }

      // method invoked by the solver
      void advop(int e)
      {
        fct_init(e); // e.g. store psi_min, psi_max in FCT

	for (int iter = 0; iter < n_iters; ++iter) 
	{
	  if (iter != 0) 
	  {
	    this->cycle(e);
            this->mem->barrier();
	    this->bcxl->fill_halos(this->mem->psi[e][this->n[e]]); // TODO: one xchng call?
	    this->bcxr->fill_halos(this->mem->psi[e][this->n[e]]);
            this->mem->barrier();

	    // calculating the antidiffusive C 
	    C_corr(iter)[0](im+h) = 
	      formulae::mpdata::antidiff(
		this->mem->psi[e][this->n[e]], 
		im, C_unco(iter)[0]
	      );

            fct_adjust_antidiff(e); // i.e. calculate C_mono=C_mono(C_corr) in FCT
          }
	  // donor-cell call
	  formulae::donorcell::op_1d(this->mem->psi[e], this->n[e], this->C(iter)[0], this->i);
	}
      }

      public:

      struct params_t {};

      // ctor
      mpdata_1d(
        typename parent_t::ctor_args_t args,
        const params_t &
      ) : 
	parent_t(args),
	im(args.i.first() - 1, args.i.last())
      {
	for (int n = 0; n < n_tmp; ++n)
          tmp[n] = &args.mem->tmp[__FILE__][n];
      }

      // memory allocation (to be called once per shared-mem node)
      static void alloc(typename parent_t::mem_t *mem, const int nx)
      {
        parent_t::alloc(mem, nx);
	for (int n = 0; n < n_tmp; ++n)
        {
	  mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>()); 
	  mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t( rng_t(0, nx-1)^h^(halo - 1) )); 
        }
      }
    };
  }; // namespace solvers
}; // namescpae libmpdataxx
