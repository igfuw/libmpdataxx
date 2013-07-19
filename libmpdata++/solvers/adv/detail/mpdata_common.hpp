/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

// TODO: template parameter:
//    enum rho_enum {rho_constant, rho_profile, rho_variable};

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <typename real_t, int n_iters_, int n_dims, int n_eqs, formulae::mpdata::opts_t opts, int minhalo>
      class mpdata_common : public detail::solver<
        real_t, n_dims, n_eqs, formulae::mpdata::n_tlev, detail::max(minhalo, formulae::mpdata::halo(opts))
      >
      {
        using parent_t = detail::solver<
          real_t, n_dims, n_eqs, formulae::mpdata::n_tlev, detail::max(minhalo, formulae::mpdata::halo(opts))
        >;

	using C_t = arrvec_t<typename parent_t::arr_t>;
 
	protected:

        // static constants
	static const int n_iters = n_iters_;
	static_assert(n_iters > 0, "n_iters <= 0");

        static const int n_tmp = n_iters > 2 ? 2 : 1; 

	// member fields
	std::array<C_t*, n_tmp> tmp;

        // methods
	C_t &C_unco(int iter)
	{   
	  return (iter == 1)  
	    ? this->mem->C 
	    : (iter % 2)  
	      ? *tmp[1]  // odd iters
	      : *tmp[0]; // even iters
	}   

	C_t &C_corr(int iter)
	{
	  return (iter  % 2)
	    ? *tmp[0]    // odd iters
	    : *tmp[1];   // even iters
	}

	virtual C_t &C(int iter)
	{
	  if (iter == 0) return this->mem->C;
	  return C_corr(iter);
	}

	// for Flux-Corrected Transport (TODO: more general names?) // TODO: move to mpdata_common
	virtual void fct_init(int e) { }
	virtual void fct_adjust_antidiff(int e, int iter) { }

        public:

        // parameters
        struct params_t {}; // TODO: abs, frac, toa, ...

        protected:

	// ctor
	mpdata_common(
	  typename parent_t::ctor_args_t args
	) : 
	  parent_t(args)
        {
	  for (int n = 0; n < n_tmp; ++n)
	    tmp[n] = &args.mem->tmp[__FILE__][n];
        }

        public:

        // 1D memory allocation
	static void alloc(typename parent_t::mem_t *mem, const int nx) 
	{   
	  parent_t::alloc(mem, nx);
	  for (int n = 0; n < n_tmp; ++n)
	    parent_t::alloc_tmp_vctr(mem, nx, __FILE__);
	}   

        // 2D version
	static void alloc(
	  typename parent_t::mem_t *mem, 
	  const int nx, const int ny
	)   
	{   
	  parent_t::alloc(mem, nx, ny);
	  for (int n = 0; n < n_tmp; ++n) 
	    parent_t::alloc_tmp_vctr(mem, nx, ny, __FILE__);
	} 

        // 3D version - TODO
      };
    }; // namespace detail

    template<typename real_t, int n_iters, int n_dims, int n_eqs, formulae::mpdata::opts_t opts, int minhalo> // TODO: reconsider arg order...
    class mpdata
    {};

    // alias names
    template <typename real_t, int n_iters, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_1d = mpdata<real_t, n_iters, 1, n_eqs, opts, minhalo>;

    template <typename real_t, int n_iters, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_2d = mpdata<real_t, n_iters, 2, n_eqs, opts, minhalo>;

    template <typename real_t, int n_iters, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_3d = mpdata<real_t, n_iters, 3, n_eqs, opts, minhalo>;
  }; // namespace solvers
}; // namescpae libmpdataxx
