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

    template<int n_iters, class mem_t, int halo = formulae::mpdata::halo>
    class mpdata_1d : public detail::solver_1d<
      mem_t, 
      formulae::mpdata::n_tlev,
      detail::max(halo, formulae::mpdata::halo)
    >
    {
      static_assert(n_iters > 0, "n_iters <= 0");

      using parent_t = detail::solver_1d< 
        mem_t, 
        formulae::mpdata::n_tlev,
        detail::max(halo, formulae::mpdata::halo)
      >;
      using arr_1d_t = typename mem_t::arr_t;

      static const int n_tmp = n_iters > 2 ? 2 : 1;

      // member fields
      arrvec_t<arr_1d_t> *tmp[n_tmp];
      rng_t im;

      protected:

      // method invoked by the solver
      void advop(int e)
      {
	for (int step = 0; step < n_iters; ++step) 
	{
	  if (step == 0) 
	    formulae::donorcell::op_1d(this->mem.psi[e], this->mem.n[e], this->mem.C[0], this->i);
	  else
	  {
	    this->cycle(e);
	    this->bcx->fill_halos(this->mem.psi[e][this->mem.n[e]]);

	    // choosing input/output for antidiff C
            const arrvec_t<arr_1d_t>
	      &C_unco = (step == 1) 
		? this->mem.C 
		: (step % 2) 
		  ? *tmp[1]  // odd steps
		  : *tmp[0], // even steps
	      &C_corr = (step  % 2) 
		? *tmp[0]    // odd steps
		: *tmp[1];   // even steps

	    // calculating the antidiffusive C 
	    C_corr[0](im+h) = 
	      formulae::mpdata::antidiff(
		this->mem.psi[e][this->mem.n[e]], 
		im, C_unco[0]
	      );

	    // donor-cell step 
	    formulae::donorcell::op_1d(this->mem.psi[e], 
	      this->mem.n[e], C_corr[0], this->i);
	  }
	}
      }

      public:

      struct params_t {};

      // ctor
      mpdata_1d(mem_t &mem, typename parent_t::bc_p &bcx, const rng_t &i, const params_t &) : 
	parent_t(mem, bcx, i), 
	im(i.first() - 1, i.last())
      {
	for (int n = 0; n < n_tmp; ++n)
          tmp[n] = &mem.tmp[std::string(__FILE__)][n];
      }

      // memory allocation (to be called once per shared-mem node)
      static void alloc(mem_t &mem, const int nx)
      {
        parent_t::alloc(mem, nx);

        const std::string file(__FILE__);

	for (int n = 0; n < n_tmp; ++n)
        {
	  mem.tmp[file].push_back(new arrvec_t<arr_1d_t>()); 
	  mem.tmp[file].back().push_back(new arr_1d_t( rng_t(0, nx-1)^h )); 
        }
      }
    };
  }; // namespace solvers
}; // namescpae advoocat
