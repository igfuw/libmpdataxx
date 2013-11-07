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

    template<typename real_t, int n_iters, int n_eqs, formulae::opts::opts_t opts, int minhalo>
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
	    this->GC_corr(iter)[0](this->im+h, this->j) = 
	      formulae::mpdata::antidiff<opts, 0>(
		this->mem->psi[e][this->n[e]], 
		this->im, this->j, this->GC_unco(iter)
	      );

	    this->GC_corr(iter)[1](this->i, this->jm+h) = 
	      formulae::mpdata::antidiff<opts, 1>(
              this->mem->psi[e][this->n[e]], 
              this->jm, this->i, this->GC_unco(iter)
	    );
 
            // filling Y halos for GC_x, and X halos for GC_y
            // TODO: document why; is it needed in the last iteration?; what about FCT?
            this->mem->barrier();
	    this->bcyl->fill_halos_sclr(this->GC_corr(iter)[0], this->i^h); // TODO: one xchng?
	    this->bcyr->fill_halos_sclr(this->GC_corr(iter)[0], this->i^h);
	    this->bcxl->fill_halos_sclr(this->GC_corr(iter)[1], this->j^h); // TODO: one xchng?
	    this->bcxr->fill_halos_sclr(this->GC_corr(iter)[1], this->j^h);
            this->mem->barrier();

            this->fct_adjust_antidiff(e, iter);

            // TODO: shouldn't the above halo-filling be repeated here?
	  }
	  // donor-cell call 
          if (!formulae::opts::isset(opts, formulae::opts::iga) || iter ==0)
	    formulae::donorcell::op_2d<opts>(
              this->mem->psi[e], 
	      this->GC(iter), 
	      this->mem->G, 
	      this->n[e], 
	      this->i, 
              this->j
            ); 
            // TODO: doing antidiff,upstream,antidiff,upstream (for each dimension separately) could help optimise memory consumption!
          else
          {
            assert(iter == 1); // infinite gauge option uses just one corrective step
	    formulae::donorcell::op_2d_iga<opts>(
              this->mem->psi[e], 
	      this->GC(iter), 
	      this->mem->G, 
	      this->n[e], 
	      this->i, 
              this->j
            ); 
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
