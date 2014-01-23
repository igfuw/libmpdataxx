/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <array>

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>  //TODO tmp

#include <libmpdata++/formulae/mpdata/formulae_mpdata_3d.hpp>
#include <libmpdata++/formulae/donorcell_formulae.hpp>
#include <libmpdata++/solvers/detail/solver_3d.hpp> // TODO: this is not used here but has to be included... tricky!
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
	typename std::enable_if<ct_params_t::n_dims == 3>::type
      > : public detail::mpdata_common<ct_params_t, minhalo>
      {
	using parent_t = detail::mpdata_common<ct_params_t, minhalo>;

	protected:

	// member fields
	rng_t im, jm, km;

	void hook_ante_loop(const int nt) 
	{   
  // TODO: same in 1D
	  parent_t::hook_ante_loop(nt);
	  if (formulae::opts::isset(ct_params_t::opts, formulae::opts::nug))
	  {
            this->bcxl->fill_halos_sclr(*this->mem->G, this->j^this->halo, this->k^this->halo); // TODO: one xchng call?
            this->bcxr->fill_halos_sclr(*this->mem->G, this->j^this->halo, this->k^this->halo);
            this->bcyl->fill_halos_sclr(*this->mem->G, this->k^this->halo, this->i^this->halo);
            this->bcyr->fill_halos_sclr(*this->mem->G, this->k^this->halo, this->i^this->halo);
            this->bczl->fill_halos_sclr(*this->mem->G, this->i^this->halo, this->j^this->halo);
            this->bczr->fill_halos_sclr(*this->mem->G, this->i^this->halo, this->j^this->halo);
	  }
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
	      this->bcxl->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->j^this->halo, this->k^this->halo); // TODO: two xchng calls? (without barriers)
	      this->bcxr->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->j^this->halo, this->k^this->halo);
	      this->bcyl->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->k^this->halo, this->i^this->halo);
	      this->bcyr->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->k^this->halo, this->i^this->halo);
	      this->bczl->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->i^this->halo, this->j^this->halo);
	      this->bczr->fill_halos_sclr(this->mem->psi[e][this->n[e]], this->i^this->halo, this->j^this->halo);
	      this->mem->barrier();

	      // calculating the antidiffusive C 
	      this->GC_corr(iter)[0](this->im+h, this->j, this->k) = 
		formulae::mpdata::antidiff<ct_params_t::opts, 0>(
                  this->mem->psi[e][this->n[e]], 
                  this->GC_unco(iter),
                  *this->mem->G,
                  this->im,
                  this->j,
                  this->k
		);

	      this->GC_corr(iter)[1](this->i, this->jm+h, this->k) = 
		formulae::mpdata::antidiff<ct_params_t::opts, 1>(
                  this->mem->psi[e][this->n[e]], 
                  this->GC_unco(iter),
                  *this->mem->G,
                  this->jm,
                  this->k,
                  this->i
	        );
	      
	      this->GC_corr(iter)[2](this->i, this->j, this->km+h) = 
		formulae::mpdata::antidiff<ct_params_t::opts, 2>(
                  this->mem->psi[e][this->n[e]], 
                  this->GC_unco(iter),
                  *this->mem->G,
                  this->km,
                  this->i,
                  this->j
	        );
   
	      // filling Y and Z halos for GC_x, X and Z halos for GC_y, X and Y
	      // halos for GC_z 
	      // TODO: document why; is it needed in the last iteration?; what about FCT?
	      this->mem->barrier();
	      this->bcyl->fill_halos_sclr(this->GC_corr(iter)[0], this->k^h, this->i^h);
	      this->bcyr->fill_halos_sclr(this->GC_corr(iter)[0], this->k^h, this->i^h);
	      this->bczl->fill_halos_sclr(this->GC_corr(iter)[0], this->i^h, this->j^h);
	      this->bczr->fill_halos_sclr(this->GC_corr(iter)[0], this->i^h, this->j^h);

	      this->bcxl->fill_halos_sclr(this->GC_corr(iter)[1], this->j^h, this->k^h);
	      this->bcxr->fill_halos_sclr(this->GC_corr(iter)[1], this->j^h, this->k^h);
	      this->bczl->fill_halos_sclr(this->GC_corr(iter)[1], this->i^h, this->j^h);
	      this->bczr->fill_halos_sclr(this->GC_corr(iter)[1], this->i^h, this->j^h);
	      
	      this->bcxl->fill_halos_sclr(this->GC_corr(iter)[2], this->j^h, this->k^h);
	      this->bcxr->fill_halos_sclr(this->GC_corr(iter)[2], this->j^h, this->k^h);
	      this->bcyl->fill_halos_sclr(this->GC_corr(iter)[2], this->k^h, this->i^h);
	      this->bcyr->fill_halos_sclr(this->GC_corr(iter)[2], this->k^h, this->i^h);
	      this->mem->barrier();

	      this->fct_adjust_antidiff(e, iter);

	      // TODO: shouldn't the above halo-filling be repeated here?
	    }
	    // donor-cell call 
	    if (!formulae::opts::isset(ct_params_t::opts, formulae::opts::iga) || iter == 0)
	    {
	      formulae::donorcell::op_3d<ct_params_t::opts>(
		this->mem->psi[e], 
		this->GC(iter), 
		*this->mem->G,
		this->n[e], 
		this->i, 
		this->j,
		this->k
	      ); 
	    }
	      // TODO: doing antidiff,upstream,antidiff,upstream (for each dimension separately) could help optimise memory consumption!
	    else
	    {
	      assert(iter == 1); // infinite gauge option uses just one corrective step
	      formulae::donorcell::op_3d_iga<ct_params_t::opts>(
		this->mem->psi[e], 
		this->GC(iter), 
		*this->mem->G,
		this->n[e], 
		this->i, 
		this->j,
		this->k
	      ); 
	    }
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
	  jm(args.j.first() - 1, args.j.last()),
	  km(args.k.first() - 1, args.k.last())
	{ }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
