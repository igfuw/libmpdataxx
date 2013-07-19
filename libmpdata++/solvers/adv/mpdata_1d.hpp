/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_1d.hpp>
#include <libmpdata++/formulae/donorcell_formulae.hpp>
#include <libmpdata++/solvers/adv/detail/solver_1d.hpp> // TODO: this is not used here but has to be included... tricky!
#include <libmpdata++/solvers/adv/detail/mpdata_common.hpp>

#include <array>

namespace libmpdataxx
{
  namespace solvers
  {
    using namespace libmpdataxx::arakawa_c;

    template<typename real_t, int n_iters, int n_eqs, formulae::mpdata::opts_t opts, int minhalo>
    class mpdata<real_t, n_iters, 1, n_eqs, opts, minhalo> : 
      public detail::mpdata_common<real_t, n_iters, 1, n_eqs, opts, minhalo>
    {
      using parent_t = detail::mpdata_common<real_t, n_iters, 1, n_eqs, opts, minhalo>;

      protected:

      rng_t im;

      // method invoked by the solver
      void advop(int e)
      {
        this->fct_init(e); // e.g. store psi_min, psi_max in FCT

	for (int iter = 0; iter < n_iters; ++iter) 
	{
	  if (iter != 0) 
	  {
	    this->cycle(e);
            this->mem->barrier();
	    this->bcxl->fill_halos_sclr(this->mem->psi[e][this->n[e]]); // TODO: one xchng call?
	    this->bcxr->fill_halos_sclr(this->mem->psi[e][this->n[e]]);
            this->mem->barrier();

	    // calculating the antidiffusive C 
	    this->C_corr(iter)[0](im+h) = 
	      formulae::mpdata::antidiff<opts>(
		this->mem->psi[e][this->n[e]], 
		im, this->C_unco(iter)[0]
	      );

            this->fct_adjust_antidiff(e, iter); // i.e. calculate C_mono=C_mono(C_corr) in FCT
          }

	  // donor-cell call
          if (!opt_set(opts, formulae::mpdata::iga) || iter == 0)
	    formulae::donorcell::op_1d(this->mem->psi[e], this->n[e], this->C(iter)[0], this->i); 
          else
          {
            assert(iter == 1); // infinite gauge option uses just one corrective step
            // TODO: move to a formulae file? 
            auto &psi = this->mem->psi[e];
            const auto &n = this->n[e];
            const auto &i = this->i;
            const auto &C = this->C(iter)[0];
            psi[n+1](i) = psi[n](i) - (C(i+h) - C(i-h)); // referred to as F(1,1,U) in the papers
          }
	}
      }

      public:

      // ctor
      mpdata(
        typename parent_t::ctor_args_t args,
        const typename parent_t::params_t &
      ) : 
	parent_t(args),
	im(args.i.first() - 1, args.i.last())
      {}
    };
  }; // namespace solvers
}; // namescpae libmpdataxx
