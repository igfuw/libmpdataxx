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
#include <libmpdata++/solvers/detail/solver_2d.hpp> // TODO: this is not used here but has to be included... tricky!
#include <libmpdata++/solvers/detail/mpdata_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      using namespace arakawa_c;

      template<typename ct_params_t, int minhalo>
      class mpdata_osc<
	ct_params_t, 
	minhalo,
	typename std::enable_if<ct_params_t::n_dims == 2>::type
      > : public detail::mpdata_common<ct_params_t, minhalo>
      {
	using parent_t = detail::mpdata_common<ct_params_t, minhalo>;

	protected:

	// member fields
	const rng_t im, jm;

	void hook_ante_loop(const int nt) 
	{   
          //  note that it's not needed for upstream
	  parent_t::hook_ante_loop(nt);
	  if (opts::isset(ct_params_t::opts, opts::nug))
	  {
	    this->bcxl->fill_halos_sclr(*this->mem->G, this->j^this->halo); // TODO: one xchng call?
	    this->bcxr->fill_halos_sclr(*this->mem->G, this->j^this->halo);
	    this->bcyl->fill_halos_sclr(*this->mem->G, this->i^this->halo);
	    this->bcyr->fill_halos_sclr(*this->mem->G, this->i^this->halo);
	  }
          
          // filling Y halos for GC_x, and X halos for GC_y
          this->bcyl->fill_halos_vctr_nrml(this->mem->GC[0], this->i^h); // TODO: one xchng?
          this->bcyr->fill_halos_vctr_nrml(this->mem->GC[0], this->i^h);
          this->bcxl->fill_halos_vctr_nrml(this->mem->GC[1], this->j^h); // TODO: one xchng?
          this->bcxr->fill_halos_vctr_nrml(this->mem->GC[1], this->j^h);
	} 

	// method invoked by the solver
	void advop(int e)
	{
	  this->fct_init(e);

	  for (int iter = 0; iter < this->n_iters; ++iter) 
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
		formulae::mpdata::antidiff<ct_params_t::opts, 0>(
		  this->mem->psi[e][this->n[e]], 
                  this->GC_unco(iter),
                  *this->mem->G,
		  this->im, 
                  this->j
		);
              assert(std::isfinite(sum(this->GC_corr(iter)[0](this->im+h, this->j))));

	      this->GC_corr(iter)[1](this->i, this->jm+h) = 
		formulae::mpdata::antidiff<ct_params_t::opts, 1>(
		  this->mem->psi[e][this->n[e]], 
		  this->GC_unco(iter),
		  *this->mem->G,
		  this->jm, 
		  this->i
		);
              assert(std::isfinite(sum(this->GC_corr(iter)[1](this->i, this->jm+h))));
   
	      // filling Y halos for GC_x, and X halos for GC_y
	      // TODO: document why; is it needed in the last iteration?; what about FCT?
	      this->mem->barrier();
	      this->bcyl->fill_halos_vctr_nrml(this->GC_corr(iter)[0], this->i^h); // TODO: one xchng?
	      this->bcyr->fill_halos_vctr_nrml(this->GC_corr(iter)[0], this->i^h);
	      this->bcxl->fill_halos_vctr_nrml(this->GC_corr(iter)[1], this->j^h); // TODO: one xchng?
	      this->bcxr->fill_halos_vctr_nrml(this->GC_corr(iter)[1], this->j^h);
	      this->mem->barrier();

	      this->fct_adjust_antidiff(e, iter);
              assert(std::isfinite(sum(this->GC_corr(iter)[0](this->im+h, this->j))));
              assert(std::isfinite(sum(this->GC_corr(iter)[1](this->i, this->jm+h))));

	      // TODO: shouldn't the above halo-filling be repeated here?
	    }

            // calculation of fluxes
            if (!opts::isset(ct_params_t::opts, opts::iga) || iter == 0)
            {
              this->flux[0](im+h, this->j) = formulae::donorcell::flux<ct_params_t::opts, 0>(
                this->mem->psi[e][this->n[e]], 
                this->GC(iter)[0], 
                im, this->j
              );
              this->flux[1](this->i, jm+h) = formulae::donorcell::flux<ct_params_t::opts, 1>(
                this->mem->psi[e][this->n[e]], 
                this->GC(iter)[1], 
                jm, this->i
              );
              this->flux_ptr = &this->flux; // TODO: if !iga this is needed only once per simulation, TODO: move to common
            }
            else
            {   
              assert(iter == 1); // infinite gauge option uses just one corrective step // TODO: not true?
              this->flux_ptr = &this->GC(iter);
            }   

            const auto &flx = (*(this->flux_ptr));

            // sanity check for input
            //assert(std::isfinite(sum(this->mem->psi[e][this->n[e]](this->ijk))));
            //assert(std::isfinite(sum(flx[0](i^h, j  ))));
            //assert(std::isfinite(sum(flx[1](i,   j^h))));

	    // donor-cell call 
	    // TODO: doing antidiff,upstream,antidiff,upstream (for each dimension separately) could help optimise memory consumption!
	    formulae::donorcell::donorcell_sum<ct_params_t::opts>(
	      this->mem->khn_tmp,
              this->ijk,
	      this->mem->psi[e][this->n[e]+1](this->ijk), 
	      this->mem->psi[e][this->n[e]  ](this->ijk), 
              flx[0](this->i+h, this->j  ),
              flx[0](this->i-h, this->j  ),
              flx[1](this->i,   this->j+h),
              flx[1](this->i,   this->j-h),
              formulae::G<ct_params_t::opts, 0>(*this->mem->G, this->i, this->j)
	    ); 

            // sanity check for output // TODO: move to common
            //assert(std::isfinite(sum(this->mem->psi[e][this->n[e]+1](this->ijk))));
	  }
	}

	public:

	// ctor
	mpdata_osc(
	  typename parent_t::ctor_args_t args,
	  const typename parent_t::rt_params_t &p
	) : 
	  parent_t(args, p),
	  im(args.i.first() - 1, args.i.last()),
	  jm(args.j.first() - 1, args.j.last())
	{ }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
