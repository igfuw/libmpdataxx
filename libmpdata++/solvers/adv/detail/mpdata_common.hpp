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
      template <typename real_t, int n_dims, formulae::mpdata::opts_t opts, int minhalo>
      class mpdata_common : public detail::solver<
        real_t, n_dims, formulae::mpdata::n_tlev, opts, detail::max(minhalo, formulae::mpdata::halo(opts))
      >
      {
        using parent_t = detail::solver<
          real_t, n_dims, formulae::mpdata::n_tlev, opts, detail::max(minhalo, formulae::mpdata::halo(opts))
        >;

	using GC_t = arrvec_t<typename parent_t::arr_t>;
 
	protected:

        // static constants
	const int n_iters;

	// member fields
	std::vector<GC_t*> tmp;

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

	// for Flux-Corrected Transport (TODO: more general names?) // TODO: move to mpdata_common
	virtual void fct_init(int e) { }
	virtual void fct_adjust_antidiff(int e, int iter) { }

        //  
        static int n_tmp(const int &n_iters)
        {
          return n_iters > 2 ? 2 : 1; 
        }

        public:

	struct params_t : parent_t::params_t
        {
          int n_iters = 2; 
        };

        protected:

	// ctor
	mpdata_common(
	  typename parent_t::ctor_args_t args,
          const params_t &p
	) : 
	  parent_t(args, p),
          n_iters(p.n_iters),
          tmp(n_tmp(n_iters))
        {
	  for (int n = 0; n < n_tmp(n_iters); ++n)
	    tmp[n] = &args.mem->tmp[__FILE__][n];

          assert(n_iters > 0);
        }

        public:

        // memory allocation
	static void alloc(
          typename parent_t::mem_t *mem, 
          const params_t &p
        ) {   
	  parent_t::alloc(mem, p);
	  for (int n = 0; n < n_tmp(p.n_iters); ++n)
	    parent_t::alloc_tmp_vctr(mem, p.span, __FILE__);
	}   
      };
    }; // namespace detail

    template<typename real_t, int n_dims, formulae::mpdata::opts_t opts, int minhalo> // TODO: reconsider arg order...
    class mpdata
    {};

    // alias names
    template <typename real_t, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_1d = mpdata<real_t, 1, opts, minhalo>;

    template <typename real_t, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_2d = mpdata<real_t, 2, opts, minhalo>;

    template <typename real_t, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_3d = mpdata<real_t, 3, opts, minhalo>;
  }; // namespace solvers
}; // namescpae libmpdataxx
