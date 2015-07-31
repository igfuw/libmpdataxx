/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <typename ct_params_t, int minhalo>
      class mpdata_common : public detail::solver<
        ct_params_t, 
        formulae::mpdata::n_tlev, 
        detail::max(minhalo, formulae::mpdata::halo(ct_params_t::opts))
      >
      {
        using parent_t = detail::solver<
          ct_params_t, 
          formulae::mpdata::n_tlev, 
          detail::max(minhalo, formulae::mpdata::halo(ct_params_t::opts))
        >;

	using GC_t = arrvec_t<typename parent_t::arr_t>;
 
	protected:

        // static constants
	const int n_iters;

	// member fields
	std::vector<GC_t*> tmp;
        GC_t &flux, *flux_ptr;

        // methods
	GC_t &GC_unco(int iter)
	{   
	  return (iter == 1)  
	    ? this->mem->GC 
	    : (iter % 2)  
	      ? *tmp[1]  // odd iters
	      : *tmp[0]; // even iters
	}   

        GC_t &GC_corr(int iter)
	{
	  return (iter  % 2)
	    ? *tmp[0]    // odd iters
	    : *tmp[1];   // even iters
	}

        virtual GC_t &GC(int iter)
	{
	  if (iter == 0) return this->mem->GC;
	  return GC_corr(iter);
	}

	// for Flux-Corrected Transport 
	virtual void fct_init(int e) { }
	virtual void fct_adjust_antidiff(int e, int iter) { }

        //  
        static int n_tmp(const int &n_iters)
        {
          return n_iters > 2 ? 2 : 1; 
        }

        public:

	struct rt_params_t : parent_t::rt_params_t
        {
          int n_iters = 2; 
        };

        protected:

	// ctor
	mpdata_common(
	  typename parent_t::ctor_args_t args,
          const rt_params_t &p
	) : 
	  parent_t(args, p),
          n_iters(p.n_iters),
          tmp(n_tmp(n_iters)),
          flux(args.mem->tmp[__FILE__][n_tmp(p.n_iters)])
        {
          assert(n_iters > 0); // TODO: throw!

	  for (int n = 0; n < n_tmp(n_iters); ++n)
	    tmp[n] = &args.mem->tmp[__FILE__][n];
        }

        public:

        // memory allocation
	static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {   
	  parent_t::alloc(mem, n_iters);
	  for (int n = 0; n < n_tmp(n_iters); ++n)
	    parent_t::alloc_tmp_vctr(mem, __FILE__);
          parent_t::alloc_tmp_vctr(mem, __FILE__); // fluxes
	}   
      };

      // partial specialisations
      template<typename ct_params_t, int minhalo, class enableif = void> 
      class mpdata_osc
      {};
    }; // namespace detail
  }; // namespace solvers
}; // namescpae libmpdataxx
