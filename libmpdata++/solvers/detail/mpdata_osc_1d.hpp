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

	rng_t im;

	void hook_ante_loop(const int nt)
	{
  // TODO: same in 2D and 3D
	  parent_t::hook_ante_loop(nt);
	  if (formulae::opts::isset(ct_params_t::opts, formulae::opts::nug))
	  {
	    this->bcxl->fill_halos_sclr(*this->mem->G); // TODO: one xchng call?
	    this->bcxr->fill_halos_sclr(*this->mem->G);
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
	      this->mem->barrier();
	      this->bcxl->fill_halos_sclr(this->mem->psi[e][this->n[e]]); // TODO: one xchng call?
	      this->bcxr->fill_halos_sclr(this->mem->psi[e][this->n[e]]);
	      this->mem->barrier();

	      // calculating the antidiffusive C 
	      this->GC_corr(iter)[0](im+h) = 
		formulae::mpdata::antidiff<ct_params_t::opts>(
		  this->mem->psi[e][this->n[e]], 
		  this->GC_unco(iter)[0],
		  *this->mem->G,
		  im
		);

	      this->fct_adjust_antidiff(e, iter); // i.e. calculate GC_mono=GC_mono(GC_corr) in FCT
	    }

	    // donor-cell call
	    if (!formulae::opts::isset(ct_params_t::opts, formulae::opts::iga) || iter == 0)
	    {
	      formulae::donorcell::op_1d<ct_params_t::opts>(
		this->mem->khn_tmp,
		this->mem->psi[e], 
		this->GC(iter)[0], 
		*this->mem->G, 
		this->n[e], 
		this->i
	      ); 
	    }
	    else
	    {
	      assert(iter == 1); // infinite gauge option uses just one corrective step
	      formulae::donorcell::op_1d_iga<ct_params_t::opts>(
		this->mem->khn_tmp,
		this->mem->psi[e], 
		this->GC(iter)[0], 
		*this->mem->G, 
		this->n[e], 
		this->i
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
	  im(args.i.first() - 1, args.i.last())
	{}
      };
    }; // namespace detail
  }; // namespace solvers
}; // namescpae libmpdataxx
