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
	const rng_t im, jm, km;
  
	void hook_ante_loop(const typename parent_t::advance_arg_t nt)
	{   
  //  note that it's not needed for upstream
	  parent_t::hook_ante_loop(nt);
	  if (opts::isset(ct_params_t::opts, opts::nug))
            this->xchng_sclr(*this->mem->G, this->i^this->halo, this->j^this->halo, this->k^this->halo);
          
          // filling Y and Z halos for GC_x, X and Z halos for GC_y, X and Y
          // halos for GC_z
          this->xchng_vctr_nrml(this->mem->GC, this->i, this->j, this->k);
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
              this->xchng(e);

	      // calculating the antidiffusive C 
	      this->GC_corr(iter)[0](this->im+h, this->j, this->k) = 
		formulae::mpdata::antidiff<ct_params_t::opts, 0>(
                  this->mem->psi[e][this->n[e]], 
                  this->GC_unco(iter),
                  this->mem->ndt_GC,
                  this->mem->ndtt_GC,
                  *this->mem->G,
                  this->im,
                  this->j,
                  this->k
		);

	      this->GC_corr(iter)[1](this->i, this->jm+h, this->k) = 
		formulae::mpdata::antidiff<ct_params_t::opts, 1>(
                  this->mem->psi[e][this->n[e]], 
                  this->GC_unco(iter),
                  this->mem->ndt_GC,
                  this->mem->ndtt_GC,
                  *this->mem->G,
                  this->jm,
                  this->k,
                  this->i
	        );
	      
	      this->GC_corr(iter)[2](this->i, this->j, this->km+h) = 
		formulae::mpdata::antidiff<ct_params_t::opts, 2>(
                  this->mem->psi[e][this->n[e]], 
                  this->GC_unco(iter),
                  this->mem->ndt_GC,
                  this->mem->ndtt_GC,
                  *this->mem->G,
                  this->km,
                  this->i,
                  this->j
	        );
              
	      // filling Y and Z halos for GC_x, X and Z halos for GC_y, X and Y halos for GC_z 
	      // needed for calculation of antidiffusive velocities in the third and subsequent
              // iterations, also needed for fct but it is done there independently hence
              // the following check
              if (!opts::isset(ct_params_t::opts, opts::fct) && iter != (this->n_iters - 1))
                this->xchng_vctr_nrml(this->GC_corr(iter), this->i, this->j, this->k);

	      this->fct_adjust_antidiff(e, iter);

	      // TODO: shouldn't the above halo-filling be repeated here?
	    }

            const auto &i(this->i), &j(this->j), &k(this->k);
            const auto &ijk(this->ijk);
            const auto &psi(this->mem->psi[e]);
            const auto &n(this->n[e]);
            auto &GC(this->GC(iter));
            using namespace formulae::donorcell;

            // calculation of fluxes
            if (!opts::isset(ct_params_t::opts, opts::iga) || iter == 0)
            {
              this->flux[0](im+h, j, k) = flux<ct_params_t::opts, 0>(psi[n], GC[0], im, j, k);
              this->flux[1](i, jm+h, k) = flux<ct_params_t::opts, 1>(psi[n], GC[1], jm, k, i);
              this->flux[2](i, j, km+h) = flux<ct_params_t::opts, 2>(psi[n], GC[2], km, i, j);
              this->flux_ptr = &this->flux; // TODO: if !iga this is needed only once per simulation, TODO: move to common
            }
            else
            {
              assert(iter == 1); // infinite gauge option uses just one corrective step // TODO: not true?
              this->flux_ptr = &GC;
            }
            
            const auto &flx = (*(this->flux_ptr));

            // sanity check for input
            assert(std::isfinite(sum(psi[n](ijk))));
            assert(std::isfinite(sum(flx[0](i^h, j,   k  ))));
            assert(std::isfinite(sum(flx[1](i,   j^h, k  ))));
            assert(std::isfinite(sum(flx[2](i,   j,   k^h))));

	    // donor-cell call 
	    // TODO: doing antidiff,upstream,antidiff,upstream (for each dimension separately) could help optimise memory consumption!
	    donorcell_sum<ct_params_t::opts>(
	      this->mem->khn_tmp,
              ijk,
	      psi[n+1](ijk), 
	      psi[n  ](ijk), 
              flx[0](i+h, j,   k  ),
              flx[0](i-h, j,   k  ),
              flx[1](i,   j+h, k  ),
              flx[1](i,   j-h, k  ),
              flx[2](i,   j,   k+h),
              flx[2](i,   j,   k-h),
              formulae::G<ct_params_t::opts, 0>(*this->mem->G, i, j, k)
	    ); 

            // sanity check for output
            assert(std::isfinite(sum(psi[n+1](ijk))));
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
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
