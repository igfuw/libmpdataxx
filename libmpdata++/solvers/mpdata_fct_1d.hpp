/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief FCT option for MPDATA as formulated in @copybrief Smolarkiewicz_and_Grabowski_1990
 */

#pragma once

#include <libmpdata++/solvers/mpdata_1d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail // TODO fct::formulae?
    {
      // TODO: document why 2
      const int fct_min_halo = 2;
    }

// TODO: could it be made once for all mpdata dimensions?
    template <typename real_t, int n_iters, int n_eqs = 1, int halo = detail::fct_min_halo> 
    class mpdata_fct_1d : public mpdata_1d<real_t, n_iters, n_eqs, detail::max(halo, detail::fct_min_halo)> 
    {
      public:

      using parent_t = mpdata_1d<real_t, n_iters, n_eqs, detail::max(halo, detail::fct_min_halo)>; 

      struct params_t : parent_t::params_t
      {
        // TODO: rho!
      };  

      private:

      typename parent_t::arr_t psi_min, psi_max;

      void fct_init(int e)
      {
// TODO: if (n_iters > 2)
std::cerr << "fct_init()" << std::endl;
        const rng_t i = this->i^1;
        const typename parent_t::arr_t psi = this->state(e);

        psi_min(i) = min(min(psi(i-1), psi(i)), psi(i+1));
        psi_max(i) = max(max(psi(i-1), psi(i)), psi(i+1));
      }

      void fct_adjust_antidiff(int e)
      {
std::cerr << "fct_adjust_antidiff()" << std::endl;
      }

      arrvec_t<typename parent_t::arr_t> &C(int iter)
      {
        //if (iter > 0) return C_mono;
        return parent_t::C(iter);
      }

      public:

      // ctor
      mpdata_fct_1d(
        typename parent_t::ctor_args_t args, 
        const params_t &p
      ) : 
        parent_t(args, p),
        psi_min(args.mem->tmp[__FILE__][0][0]),
        psi_max(args.mem->tmp[__FILE__][0][1])
      {}   

      // 1D version
      static void alloc(typename parent_t::mem_t *mem, const int nx)
      {
	parent_t::alloc(mem, nx);

        const rng_t i(0, nx-1);

	mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
	for (int e = 0; e < 2; ++e)
	  mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t(i^parent_t::halo));  // TODO: halo=1 would be enough!
      }
    };
  }; // namespace solvers
}; // namespace libmpdataxx
