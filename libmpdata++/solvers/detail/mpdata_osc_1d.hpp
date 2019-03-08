/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_1d.hpp>
#include <libmpdata++/formulae/donorcell_formulae.hpp>
#include <libmpdata++/solvers/detail/solver_1d.hpp> // TODO: this is not used here but has to be included... tricky!
#include <libmpdata++/solvers/detail/mpdata_common.hpp>

#include <array>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      using namespace libmpdataxx::arakawa_c;

      template<typename ct_params_t, int minhalo>
      class mpdata_osc<
	ct_params_t, 
	minhalo,
	typename std::enable_if<ct_params_t::n_dims == 1>::type
      > : public detail::mpdata_common<ct_params_t, minhalo>
      {
	using parent_t = detail::mpdata_common<ct_params_t, minhalo>;

	protected:

	const rng_t im;

	void hook_ante_loop(const typename parent_t::advance_arg_t nt)
	{
  //  note that it's not needed for upstream
	  parent_t::hook_ante_loop(nt);
	  if (opts::isset(ct_params_t::opts, opts::nug))
            this->xchng_sclr(*this->mem->G); 

          // set time derivatives of GC to zero
          // needed for stationary flows prescribed using the advector method
          if (opts::isset(ct_params_t::opts, opts::div_3rd_dt) || opts::isset(ct_params_t::opts, opts::div_3rd))
          {
            this->mem->ndt_GC[0](this->im + h) = 0;
            
            this->mem->ndtt_GC[0](this->im + h) = 0;

            this->xchng_vctr_alng(this->mem->ndt_GC);
            this->xchng_vctr_alng(this->mem->ndtt_GC);
          }
	}

	// method invoked by the solver
	void advop(int e)
	{
	  this->fct_init(e); // e.g. store psi_min, psi_max in FCT

	  for (int iter = 0; iter < this->n_iters; ++iter) 
	  {
	    if (iter != 0) 
	    {
	      this->cycle(e); // cycles subdomain's "n", and global "n" if it's the last equation
              this->xchng(e);

	      // calculating the antidiffusive C 
              formulae::mpdata::antidiff<ct_params_t::opts,
                                         static_cast<sptl_intrp_t>(ct_params_t::sptl_intrp),
                                         static_cast<tmprl_extrp_t>(ct_params_t::tmprl_extrp)>(
                this->GC_corr(iter)[0],
                this->mem->psi[e][this->n[e]], 
                this->GC_unco(iter),
                this->mem->ndt_GC,
                this->mem->ndtt_GC,
                *this->mem->G,
                im
              );

              // needed with the dfl option
              // if we aren't in the last iteration and fct is not set
              if (opts::isset(ct_params_t::opts, opts::dfl) &&
                  iter != (this->n_iters - 1) &&
                  !opts::isset(ct_params_t::opts, opts::fct))
              {
                this->xchng_vctr_alng(this->GC_corr(iter));
              }

	      this->fct_adjust_antidiff(e, iter); // i.e. calculate GC_mono=GC_mono(GC_corr) in FCT
	    }

	    // calculation of fluxes
	    if (!opts::isset(ct_params_t::opts, opts::iga) || iter == 0)
	    {
              this->flux[0](im+h) = formulae::donorcell::make_flux<ct_params_t::opts>(
                this->mem->psi[e][this->n[e]],
                this->GC(iter)[0], 
                im
              );
              this->flux_ptr = &this->flux; // TODO: if !iga this is needed only once per simulation, TODO: move to common
	    }
	    else
	    {
	      assert(iter == 1); // infinite gauge option uses just one corrective step // TODO: not true?
              this->flux_ptr = &this->GC(iter); // TODO: move to common
	    }

            // sanity checks for input // TODO: move to common
            //assert(std::isfinite(sum(psi[this->n[e]](this->ijk)))); 
            //assert(std::isfinite(sum(flux_ref[0](i^h))));

	    // donor-cell call // TODO: could be made common for 1D/2D/3D
            formulae::donorcell::donorcell_sum<ct_params_t::opts>(
              this->mem->khn_tmp,
              this->ijk,
              this->mem->psi[e][this->n[e]+1](this->ijk),
              this->mem->psi[e][this->n[e]  ](this->ijk),
              (*(this->flux_ptr))[0](this->i+h),
              (*(this->flux_ptr))[0](this->i-h),
              formulae::G<ct_params_t::opts>(*this->mem->G, this->i)
            );

            // sanity checks for output // TODO: move to common
            //assert(std::isfinite(sum(psi[this->n[e]+1](this->ijk)))); 
	  }
	}

        // performs advection of a given field using the donorcell scheme
        // and stores the result in the same field
        // useful for advecting right-hand-sides etc
        void self_advec_donorcell(typename parent_t::arr_t &field)
        {
          const auto &i(this->i);
          auto &GC(this->mem->GC);
          using namespace formulae::donorcell;
 
          this->xchng_sclr(field, this->ijk);
 
          // calculation of fluxes
          this->flux[0](im+h) = make_flux<ct_params_t::opts>(field, GC[0], im);
 
          // sanity check for input
          assert(std::isfinite(sum(field(i))));
          assert(std::isfinite(sum(this->flux[0](i^h))));
 
          // donor-cell call 
          donorcell_sum<ct_params_t::opts>(
            this->mem->khn_tmp,
            i,
            field(i),
            field(i),
            this->flux[0](i+h),
            this->flux[0](i-h),
            formulae::G<ct_params_t::opts>(*this->mem->G, i)
          );

          // sanity check for output
          assert(std::isfinite(sum(field(i))));
        }


	public:

	// ctor
	mpdata_osc(
	  typename parent_t::ctor_args_t args,
	  const typename parent_t::rt_params_t &p
	) : 
	  parent_t(args, p),
	  im(args.i.first() - 1, args.i.last())
	{}
      };
    } // namespace detail
  } // namespace solvers
} // namescpae libmpdataxx
