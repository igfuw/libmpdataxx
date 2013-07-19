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

      // TODO: document why 2
      const int fct_min_halo = 2; // TODO move to fct::formulae?

      template <typename real_t, int n_iters, int n_dims, int n_eqs, formulae::mpdata::opts_t opts, int minhalo>
      class mpdata_fct_common : public mpdata<
        real_t, n_iters, n_dims, n_eqs, opts, detail::max(minhalo, fct_min_halo)
      >
      {
        using parent_t = mpdata<
          real_t, n_iters, n_dims, n_eqs, opts, detail::max(minhalo, fct_min_halo)
        >;

        static_assert(parent_t::n_iters > 1, "FCT is defined for MPDATA with a corrective iteration (not for donorcell)");
 
        protected:

        // member fields
	typename parent_t::arr_t psi_min, psi_max; 
	arrvec_t<typename parent_t::arr_t> C_mono; 

	arrvec_t<typename parent_t::arr_t> &C(int iter) 
	{
	  if (iter > 0) return C_mono;
	  return parent_t::C(iter);
	}

        public:

	struct params_t : parent_t::params_t
	{
	  // TODO: rho!
	};

        protected:

        // ctor
	mpdata_fct_common(
	  typename parent_t::ctor_args_t args,
	  const params_t &p
	) : 
          parent_t(args, p),
	  psi_min(args.mem->tmp[__FILE__][0][0]),
	  psi_max(args.mem->tmp[__FILE__][0][1]),
	   C_mono(args.mem->tmp[__FILE__][1])
        {}

        public:

	// 1D version
	static void alloc(typename parent_t::mem_t *mem, const int nx)
	{
	  parent_t::alloc(mem, nx);
	  parent_t::alloc_tmp_sclr(mem, nx, __FILE__, 2); // psi_min and psi_max
	  parent_t::alloc_tmp_vctr(mem, nx, __FILE__);    // C_mono
	}

        // 2D version 
	static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
	{
	  parent_t::alloc(mem, nx, ny);
	  parent_t::alloc_tmp_sclr(mem, nx, ny, __FILE__, 2); // psi_min and psi_max
	  parent_t::alloc_tmp_vctr(mem, nx, ny, __FILE__);    // C_mono
	}

        // 3D version - TODO
      };
    }; // namespace detail

    template<typename real_t, int n_iters, int n_dims, int n_eqs, formulae::mpdata::opts_t opts, int minhalo> // TODO: reconsider arg order...
    class mpdata_fct
    {};

    // alias names
    template <typename real_t, int n_iters, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_fct_1d = mpdata_fct<real_t, n_iters, 1, n_eqs, opts, minhalo>;

    template <typename real_t, int n_iters, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_fct_2d = mpdata_fct<real_t, n_iters, 2, n_eqs, opts, minhalo>;

    template <typename real_t, int n_iters, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_fct_3d = mpdata_fct<real_t, n_iters, 3, n_eqs, opts, minhalo>;
  }; // namespace solvers
}; // namescpae libmpdataxx
