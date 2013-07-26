/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <array>

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>  //TODO tmp

#include <libmpdata++/formulae/mpdata/formulae_mpdata_2d.hpp>
#include <libmpdata++/formulae/donorcell_formulae.hpp>
#include <libmpdata++/solvers/adv/detail/solver_2d.hpp> // TODO: this is not used here but has to be included... tricky!
#include <libmpdata++/solvers/adv/detail/mpdata_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    using namespace arakawa_c;

    template<typename real_t, int n_iters, int n_eqs, formulae::mpdata::opts_t opts, int minhalo>
    class mpdata<real_t, n_iters, 2, n_eqs, opts, minhalo> : 
      public detail::mpdata_common<real_t, n_iters, 2, n_eqs, opts, minhalo>
    {
      using parent_t = detail::mpdata_common<real_t, n_iters, 2, n_eqs, opts, minhalo>;

      protected:

      // member fields
      rng_t im, jm;

      // method invoked by the solver
      void advop(int e)
      {
        this->fct_init(e);

	for (int iter = 0; iter < n_iters; ++iter) 
	{
          if (iter != 0)
	  {
	    this->cycle(e);
            this->mem->barrier();
	    this->bcxl->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->j^this->halo); // TODO: two xchng calls? (without barriers)
	    this->bcxr->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->j^this->halo);
	    this->bcyl->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->i^this->halo);
	    this->bcyr->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->i^this->halo);
            this->mem->barrier();

	    // calculating the antidiffusive C 
	    this->C_corr(iter)[0](this->im+h, this->j) = 
	      formulae::mpdata::antidiff<opts, 0>(
		this->mem->psi[e][this->n[e]], 
		this->im, this->j, this->C_unco(iter)
	      );

	    this->C_corr(iter)[1](this->i, this->jm+h) = 
	      formulae::mpdata::antidiff<opts, 1>(
              this->mem->psi[e][this->n[e]], 
              this->jm, this->i, this->C_unco(iter)
	    );
 
            // filling Y halos for C_x, and X halos for C_y
            // TODO: document why; is it needed in the last iteration?; what about FCT?
            this->mem->barrier();
	    this->bcyl->fill_halos_sclr(this->C_corr(iter)[0], this->i^h); // TODO: one xchng?
	    this->bcyr->fill_halos_sclr(this->C_corr(iter)[0], this->i^h);
	    this->bcxl->fill_halos_sclr(this->C_corr(iter)[1], this->j^h); // TODO: one xchng?
	    this->bcxr->fill_halos_sclr(this->C_corr(iter)[1], this->j^h);
            this->mem->barrier();

            this->fct_adjust_antidiff(e, iter);

            // TODO: shouldn't the above halo-filling be repeated here?
	  }
	  // donor-cell call 
          if (!opt_set(opts, formulae::mpdata::iga) || iter ==0)
	    formulae::donorcell::op_2d(this->mem->psi[e], this->n[e], this->C(iter), this->i, this->j); 
            // TODO: doing antidiff,upstream,antidiff,upstream (for each dimension separately) could help optimise memory consumption!
          else
          {
            assert(iter == 1); // infinite gauge option uses just one corrective step
            // TODO: move to a formulae file? 
            auto &psi = this->mem->psi[e];
            const auto &n = this->n[e];
            const auto &i = this->i;
            const auto &j = this->j;
            const auto &C = this->C(iter);
            psi[n+1](i, j) = psi[n](i, j) - (
              (C[0](i+h, j) - C[0](i-h, j)) +
              (C[1](i, j+h) - C[1](i, j-h))
            ); // referred to as F(1,1,U) in the papers
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
	im(args.i.first() - 1, args.i.last()),
	jm(args.j.first() - 1, args.j.last())
      { }
    };
  }; // namespace solvers
}; // namespace libmpdataxx
