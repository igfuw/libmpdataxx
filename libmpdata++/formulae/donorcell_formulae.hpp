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
      ) return_macro(, 
        pospart<opts>(GC) * psi_l +
        negpart<opts>(GC) * psi_r
      ) 

      template <opts_t opts, class arr_1d_t>
      inline auto flux_rght( 
	const arr_1d_t &psi, 
	const arr_1d_t &GC, 
	const rng_t &i
      ) return_macro(,
	-F<opts>(
	  psi(i  ), 
	  psi(i+1), 
	   GC(i+h)
	) 
      )

      template <opts_t opts, class arr_1d_t>
      inline auto flux_left( 
	const arr_1d_t &psi, 
	const arr_1d_t &GC, 
	const rng_t &i
      ) return_macro(,
	F<opts>(
	  psi(i-1), 
	  psi(i  ), 
	   GC(i-h)
	)
      )

      template<opts_t opts, int d, class arr_2d_t>  
      inline auto flux_rght( 
	const arr_2d_t &psi, 
	const arr_2d_t &GC, 
	const rng_t &i, 
	const rng_t &j
      ) return_macro(,
	- F<opts>(
	  psi(pi<d>(i,   j)), 
	  psi(pi<d>(i+1, j)), 
	   GC(pi<d>(i+h, j))
	)
      )

      template<opts_t opts, int d, class arr_2d_t>  
      inline auto flux_left( 
	const arr_2d_t &psi, 
	const arr_2d_t &GC, 
	const rng_t &i, 
	const rng_t &j
      ) return_macro(,
	F<opts>(
	  psi(pi<d>(i-1, j)), 
	  psi(pi<d>(i,   j)), 
	   GC(pi<d>(i-h, j))
	)
      )

      template<opts_t opts, int d, class arr_3d_t>  
      inline auto flux_rght( 
	const arr_3d_t &psi, 
	const arr_3d_t &GC, 
	const rng_t &i, 
	const rng_t &j,
	const rng_t &k
      ) return_macro(,
	- F<opts>(
	  psi(pi<d>(i,   j, k)), 
	  psi(pi<d>(i+1, j, k)), 
	   GC(pi<d>(i+h, j, k))
	) 
      )

      template<opts_t opts, int d, class arr_3d_t>  
      inline auto flux_left( 
	const arr_3d_t &psi, 
	const arr_3d_t &GC, 
	const rng_t &i, 
	const rng_t &j,
	const rng_t &k
      ) return_macro(,
	F<opts>(
	  psi(pi<d>(i-1, j, k)), 
	  psi(pi<d>(i,   j, k)), 
	   GC(pi<d>(i-h, j, k))
	)
      )

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
	  psi_new = psi_old + (flx_1 + flx_2) / g;
        }
        else
        {
          kahan_zro(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new);
          kahan_add(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new, psi_old);
          kahan_add(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new, flx_1 / g);
          kahan_add(khn_tmp[0](i), khn_tmp[1](i), khn_tmp[2](i), psi_new, flx_2 / g);
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
	  psi_new = psi_old + ((flx_1 + flx_2) + (flx_3 + flx_4)) / g;
        }
        else
        { 
          kahan_zro(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, psi_old);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, flx_1 / g);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, flx_2 / g);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, flx_3 / g);
          kahan_add(khn_tmp[0](ij), khn_tmp[1](ij), khn_tmp[2](ij), psi_new, flx_4 / g);
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
	  psi_new = psi_old + ((flx_1 + flx_2) + (flx_3 + flx_4) + (flx_5 + flx_6)) / g;
        } 
        else
        {
          kahan_zro(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, psi_old);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, flx_1 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, flx_2 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, flx_3 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, flx_4 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, flx_5 / g);
          kahan_add(khn_tmp[0](ijk), khn_tmp[1](ijk), khn_tmp[2](ijk), psi_new, flx_6 / g);
        }
      }

      template <opts_t opts, class arr_1d_t>
      inline void op_1d(
	const arrvec_t<arr_1d_t> &khn_tmp,
	const arrvec_t<arr_1d_t> &psi, 
	const arr_1d_t &GC, 
	const arr_1d_t &G, 
	const int n,
	const rng_t &i
      ) { 
        assert(std::isfinite(sum(psi[n](i))));
        assert(std::isfinite(sum(GC(i^h))));

        donorcell_sum<opts>(
          khn_tmp,
          idx_t<1>(i),
	  psi[n+1](i), 
          psi[n](i),
	  flux_left<opts>(psi[n], GC, i),
	  flux_rght<opts>(psi[n], GC, i),
          formulae::G<opts>(G, i)
        );

        assert(std::isfinite(sum(psi[n+1](i))));
      }

      // infinite-gauge version (referred to as F(1,1,U) in the papers)
      template <opts_t opts, class arr_1d_t>
      inline void op_1d_iga(
	const arrvec_t<arr_1d_t> &khn_tmp,
	const arrvec_t<arr_1d_t> &psi, 
	const arr_1d_t &GC, 
	const arr_1d_t &G, 
	const int n,
	const rng_t &i
      ) { 
        assert(std::isfinite(sum(psi[n](i))));
        assert(std::isfinite(sum(GC(i^h))));

        donorcell_sum<opts>(
          khn_tmp,
          idx_t<1>(i),
	  psi[n+1](i),
          psi[n](i),
          -GC(i+h),
           GC(i-h),
          formulae::G<opts>(G, i)
        );

        assert(std::isfinite(sum(psi[n+1](i))));
      }

      template <opts_t opts, class arr_2d_t>
      inline void op_2d(
	const arrvec_t<arr_2d_t> &khn_tmp,
	const arrvec_t<arr_2d_t> &psi,
	const arrvec_t<arr_2d_t> &GC, 
	const arr_2d_t &G, 
        const int n,
	const rng_t &i, const rng_t &j
      ) { 
        assert(std::isfinite(sum(psi[n](i,   j  ))));
        assert(std::isfinite(sum(GC[0 ](i^h, j  ))));
        assert(std::isfinite(sum(GC[1 ](i,   j^h))));

        donorcell_sum<opts>(
          khn_tmp, 
          idx_t<2>({i, j}),
	  psi[n+1](i,j),
          psi[n](i,j),
	  flux_left<opts, 0>(psi[n], GC[0], i, j), 
	  flux_rght<opts, 0>(psi[n], GC[0], i, j),
	  flux_left<opts, 1>(psi[n], GC[1], j, i),
	  flux_rght<opts, 1>(psi[n], GC[1], j, i),          
          formulae::G<opts, 0>(G, i, j)
        );

        assert(std::isfinite(sum(psi[n+1](i,j))));
      }

      // infinite-gauge version (referred to as F(1,1,U) in the papers)
      template <opts_t opts, class arr_2d_t>
      inline void op_2d_iga(
	const arrvec_t<arr_2d_t> &khn_tmp,
	const arrvec_t<arr_2d_t> &psi,
	const arrvec_t<arr_2d_t> &GC, 
	const arr_2d_t &G, 
        const int n,
	const rng_t &i, const rng_t &j
      ) { 
        assert(std::isfinite(sum(psi[n](i,   j  ))));
        assert(std::isfinite(sum(GC[0 ](i^h, j  ))));
        assert(std::isfinite(sum(GC[1 ](i,   j^h))));

        donorcell_sum<opts>(
          khn_tmp,
          idx_t<2>({i, j}),
	  psi[n+1](i,j),
          psi[n](i,j),
          -GC[0](i+h, j),
           GC[0](i-h, j),
          -GC[1](i, j+h),
           GC[1](i, j-h),
          formulae::G<opts, 0>(G, i, j)
        );

        assert(std::isfinite(sum(psi[n+1](i,j))));
      }

      template <opts_t opts, class arr_3d_t>
      inline void op_3d(
	const arrvec_t<arr_3d_t> &khn_tmp,
	const arrvec_t<arr_3d_t> &psi, 
	const arrvec_t<arr_3d_t> &GC, 
	const arr_3d_t &G, 
        const int n,
	const rng_t &i, const rng_t &j, const rng_t &k
      ) { 
        assert(std::isfinite(sum(psi[n](i,   j,   k  ))));
        assert(std::isfinite(sum(GC[0 ](i^h, j  , k  ))));
        assert(std::isfinite(sum(GC[1 ](i,   j^h, k  ))));
        assert(std::isfinite(sum(GC[2 ](i,   j,   k^h))));

        donorcell_sum<opts>(
          khn_tmp,
          idx_t<3>({i,j,k}),
	  psi[n+1](i,j,k),
          psi[n](i,j,k),
	  flux_left<opts, 0>(psi[n], GC[0], i, j, k),
	  flux_rght<opts, 0>(psi[n], GC[0], i, j, k),
	  flux_left<opts, 1>(psi[n], GC[1], j, k, i),
	  flux_rght<opts, 1>(psi[n], GC[1], j, k, i),
	  flux_left<opts, 2>(psi[n], GC[2], k, i, j),
	  flux_rght<opts, 2>(psi[n], GC[2], k, i, j),
          formulae::G<opts, 0>(G, i, j, k)
        );

        assert(std::isfinite(sum(psi[n+1](i,j,k))));
      }
      
      // infinite-gauge version (referred to as F(1,1,U) in the papers)
      template <opts_t opts, class arr_3d_t>
      inline void op_3d_iga(
	const arrvec_t<arr_3d_t> &khn_tmp,
	const arrvec_t<arr_3d_t> &psi,
	const arrvec_t<arr_3d_t> &GC, 
	const arr_3d_t &G, 
        const int n,
	const rng_t &i, const rng_t &j, const rng_t &k
      ) { 
        assert(std::isfinite(sum(psi[n](i,   j,   k  ))));
        assert(std::isfinite(sum(GC[0 ](i^h, j  , k  ))));
        assert(std::isfinite(sum(GC[1 ](i,   j^h, k  ))));
        assert(std::isfinite(sum(GC[2 ](i,   j,   k^h))));

        donorcell_sum<opts>(
          khn_tmp,
          idx_t<3>({i,j,k}),
	  psi[n+1](i,j,k),
          psi[n](i,j,k),
	  -GC[0](i+h, j, k), 
	   GC[0](i-h, j, k),
	  -GC[1](i, j+h, k),
	   GC[1](i, j-h, k),
	  -GC[2](i, j, k+h),
	   GC[2](i, j, k-h),
          formulae::G<opts, 0>(G, i, j, k)
        ); 

        assert(std::isfinite(sum(psi[n+1](i,j,k))));
      }

    }; // namespace donorcell 
  }; // namespace formulae
}; // namespace libmpdataxx
