/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/idxperm.hpp>
#include <libmpdata++/formulae/common.hpp>
#include <libmpdata++/formulae/kahan_sum.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace donorcell
    {
      using namespace arakawa_c;
      using idxperm::pi;
      using opts::opts_t;

      const int n_tlev = 2, halo = 1;

      template<opts_t opts, class T1, class T2, class T3> 
      inline auto F(
	const T1 &psi_l, const T2 &psi_r, const T3 &GC
      )
      {
	return return_helper<rng_t>(
          pospart<opts BOOST_PP_COMMA() rng_t>(GC) * psi_l +
          negpart<opts BOOST_PP_COMMA() rng_t>(GC) * psi_r
        );
      }

      template <opts_t opts, class arr_1d_t>
      inline auto make_flux( 
	const arr_1d_t &psi, 
	const arr_1d_t &GC, 
	const rng_t &i
      )
      {
	return return_helper<rng_t>(F<opts>(
	  psi(i  ), 
	  psi(i+1), 
	  GC(i+h)
        ));
      }

      template<opts_t opts, int d, class arr_2d_t>  
      inline auto make_flux( 
	const arr_2d_t &psi, 
	const arr_2d_t &GC, 
	const rng_t &i, 
	const rng_t &j
      ) 
      {
	return return_helper<rng_t>(F<opts>(
	  psi(pi<d>(i,   j)), 
	  psi(pi<d>(i+1, j)), 
	   GC(pi<d>(i+h, j))
	));
      }

      template<opts_t opts, int d, class arr_3d_t>  
      inline auto make_flux( 
	const arr_3d_t &psi, 
	const arr_3d_t &GC, 
	const rng_t &i, 
	const rng_t &j,
	const rng_t &k
      )
      {
	return return_helper<rng_t>(F<opts>(
	  psi(pi<d>(i,   j, k)), 
	  psi(pi<d>(i+1, j, k)), 
	   GC(pi<d>(i+h, j, k))
	));
      }

      template <opts_t opts, class a_t, class f1_t, class f2_t, class g_t>
      inline void donorcell_sum(
        const arrvec_t<a_t> &khn_tmp,
        const idx_t<1> i,
        a_t psi_new, 
        const a_t &psi_old,
        const f1_t &flx_1,
        const f2_t &flx_2,
        const g_t &g
      )
      {
        if (!opts::isset(opts, opts::khn))
        {
	  psi_new = psi_old + (-flx_1 + flx_2) / g;
        }
        else
        {
          kahan_zro(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new);
          kahan_add(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new, psi_old);
          kahan_add(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new, -flx_1 / g);
          kahan_add(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new,  flx_2 / g);
        }
      }

      template <opts_t opts, class a_t, class f1_t, class f2_t, class f3_t, class f4_t, class g_t>
      inline void donorcell_sum(
        const arrvec_t<a_t> &khn_tmp,
        const idx_t<2> ij,
        a_t psi_new,
        const a_t &psi_old,
        const f1_t &flx_1,
        const f2_t &flx_2,
        const f3_t &flx_3,
        const f4_t &flx_4,
        const g_t &g
      )
      {
        if (!opts::isset(opts, opts::khn))
        {
	  // note: the parentheses are intended to minimise chances of numerical errors
	  psi_new = psi_old + ((-flx_1 + flx_2) + (-flx_3 + flx_4)) / g;
        }
        else
        { 
          kahan_zro(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, psi_old);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, -flx_1 / g);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new,  flx_2 / g);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, -flx_3 / g);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new,  flx_4 / g);
        }
      }

      template <opts_t opts, class a_t, class f1_t, class f2_t, class f3_t, class f4_t, class f5_t, class f6_t, class g_t>
      inline void donorcell_sum(
        const arrvec_t<a_t> &khn_tmp,
        const idx_t<3> ijk,
        a_t psi_new,
        const a_t &psi_old,
        const f1_t &flx_1,
        const f2_t &flx_2,
        const f3_t &flx_3,
        const f4_t &flx_4,
        const f5_t &flx_5,
        const f6_t &flx_6,
        const g_t &g
      )
      {
        if (!opts::isset(opts, opts::khn))
        {
	  // note: the parentheses are intended to minimise chances of numerical errors
	  psi_new = psi_old + ((-flx_1 + flx_2) + (-flx_3 + flx_4) + (-flx_5 + flx_6)) / g;
        } 
        else
        {
          kahan_zro(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, psi_old);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, -flx_1 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new,  flx_2 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, -flx_3 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new,  flx_4 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, -flx_5 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new,  flx_6 / g);
        }
      }

    } // namespace donorcell 
  } // namespace formulae
} // namespace libmpdataxx
