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

      const int fct_min_halo = 2; // TODO move to fct::formulae? & document why 2

      template <typename real_t, int n_dims, int n_eqs, formulae::mpdata::opts_t opts, int minhalo>
      class mpdata_fct_common : public mpdata<
        real_t, n_dims, n_eqs, opts, detail::max(minhalo, fct_min_halo)
      >
      {
        using parent_t = mpdata<
          real_t, n_dims, n_eqs, opts, detail::max(minhalo, fct_min_halo)
        >;

        protected:

        // member fields
	typename parent_t::arr_t psi_min, psi_max; 
	arrvec_t<typename parent_t::arr_t> GC_mono; 

	arrvec_t<typename parent_t::arr_t> &GC(int iter) 
	{
	  if (iter > 0) return GC_mono;
	  return parent_t::GC(iter);
	}

        protected:

        // ctor
	mpdata_fct_common(
	  typename parent_t::ctor_args_t args,
	  const typename parent_t::params_t &p
	) : 
          parent_t(args, p),
	  psi_min(args.mem->tmp[__FILE__][0][0]),
	  psi_max(args.mem->tmp[__FILE__][0][1]),
	  GC_mono(args.mem->tmp[__FILE__][1])
        {
          assert(parent_t::n_iters > 1 && "FCT is defined for MPDATA with a corrective iteration (not for donorcell)");
        }

        public:

	static void alloc(typename parent_t::mem_t *mem, const typename parent_t::params_t &p)
	{
	  parent_t::alloc(mem, p);
	  parent_t::alloc_tmp_sclr(mem, p.span, __FILE__, 2); // psi_min and psi_max
	  parent_t::alloc_tmp_vctr(mem, p.span, __FILE__);    // GC_mono
	}
      };
    }; // namespace detail

    template<typename real_t, int n_dims, int n_eqs, formulae::mpdata::opts_t opts, int minhalo> // TODO: reconsider arg order...
    class mpdata_fct
    {};

    // alias names
    template <typename real_t, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_fct_1d = mpdata_fct<real_t, 1, n_eqs, opts, minhalo>;

    template <typename real_t, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_fct_2d = mpdata_fct<real_t, 2, n_eqs, opts, minhalo>;

    template <typename real_t, int n_eqs = 1, formulae::mpdata::opts_t opts = 0, int minhalo = 0>
    using mpdata_fct_3d = mpdata_fct<real_t, 3, n_eqs, opts, minhalo>;
  }; // namespace solvers
}; // namescpae libmpdataxx
