/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "../formulae/mpdata/formulae_mpdata_1d.hpp"
#include "../formulae/donorcell_formulae.hpp"
#include "detail/solver_1d.hpp"
#include <unordered_map>

// TODO: an mpdata_common class?

namespace advoocat
{
  namespace solvers
  {
    using namespace advoocat::arakawa_c;

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

      using parent_t = detail::solver_1d< 
        real_t,
        1,
        n_eqs,
        formulae::mpdata::n_tlev,
        detail::max(halo, formulae::mpdata::halo)
      >;

      static const int n_tmp = n_iters > 2 ? 2 : 1;

      // member fields
      arrvec_t<typename parent_t::arr_t> *tmp[n_tmp];
      rng_t im;

      protected:

      // method invoked by the solver
      void advop(int e)
      {
	for (int step = 0; step < n_iters; ++step) 
	{
	  if (step == 0) 
	    formulae::donorcell::op_1d(this->mem->psi[e], this->n[e], this->mem->C[0], this->i);
	  else
	  {
	    this->cycle(e);
            this->mem->barrier();
	    this->bcxl->fill_halos(this->mem->psi[e][this->n[e]]); // TODO: one xchng call?
	    this->bcxr->fill_halos(this->mem->psi[e][this->n[e]]);
            this->mem->barrier();

	    // choosing input/output for antidiff C
            const arrvec_t<typename parent_t::arr_t>
	      &C_unco = (step == 1) 
		? this->mem->C 
		: (step % 2) 
		  ? *tmp[1]  // odd steps
		  : *tmp[0], // even steps
	      &C_corr = (step  % 2) 
		? *tmp[0]    // odd steps
		: *tmp[1];   // even steps

	    // calculating the antidiffusive C 
	    C_corr[0](im+h) = 
	      formulae::mpdata::antidiff(
		this->mem->psi[e][this->n[e]], 
		im, C_unco[0]
	      );

	    // donor-cell step 
	    formulae::donorcell::op_1d(this->mem->psi[e], 
	      this->n[e], C_corr[0], this->i);
	  }
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
	  mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t( rng_t(0, nx-1)^h )); 
        }
      }
    };
  }; // namespace solvers
}; // namescpae advoocat
