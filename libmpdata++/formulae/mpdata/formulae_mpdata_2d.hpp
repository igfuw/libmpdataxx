/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {

      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(  // positive sign scalar version
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && opt_set(opts, pds)>::type* = 0
      ) return_macro(,
	frac<opts>(
	  psi(pi<d>(i+1, j)) - psi(pi<d>(i, j))
	  ,// --------------------------------
	  psi(pi<d>(i+1, j)) + psi(pi<d>(i, j))
	)
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(  // variable-sign scalar version
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0
      ) return_macro(,
	frac<opts>(
	  abs(psi(pi<d>(i+1, j))) - abs(psi(pi<d>(i, j)))
	  ,// -------------------------------------------
	  abs(psi(pi<d>(i+1, j))) + abs(psi(pi<d>(i, j)))
	) 
      ) 

      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(  // inf. gauge option
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<opt_set(opts, iga)>::type* = 0
      ) return_macro(,
        (
	  psi(pi<d>(i+1, j)) - psi(pi<d>(i, j))
        ) / ( //---------------------
	  1 + 1
        )
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B( // positive sign signal
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && opt_set(opts, pds)>::type* = 0
      ) return_macro(,
	frac<opts>( 
	  psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) - psi(pi<d>(i+1, j-1)) - psi(pi<d>(i, j-1))
	  ,// --------------------------------------------------------------------------------
	  psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) + psi(pi<d>(i+1, j-1)) + psi(pi<d>(i, j-1))
	) / 2
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B( // variable-sign signal
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0
      ) return_macro(,
	frac<opts>( 
	  abs(psi(pi<d>(i+1, j+1))) + abs(psi(pi<d>(i, j+1))) - abs(psi(pi<d>(i+1, j-1))) - abs(psi(pi<d>(i, j-1)))
	  ,// ----------------------------------------------------------------------------------------------------
	  abs(psi(pi<d>(i+1, j+1))) + abs(psi(pi<d>(i, j+1))) + abs(psi(pi<d>(i+1, j-1))) + abs(psi(pi<d>(i, j-1)))
	) / 2
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B( // inf. gauge
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<opt_set(opts, iga)>::type* = 0
      ) return_macro(,
        (
	  psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) - psi(pi<d>(i+1, j-1)) - psi(pi<d>(i, j-1))
	) / (  // --------------------------------------------------------------------------------
	  1 + 1 + 1 +1
        ) / 2
      )

      template<int d, class arr_2d_t>
      inline auto C_bar( // caution proper call looks like C_bar<dim>(C[dim-1], i, j) - note dim vs dim-1
        const arr_2d_t &C, 
        const rng_t &i, 
        const rng_t &j
      ) return_macro(,
        (
          C(pi<d>(i+1, j+h)) + 
          C(pi<d>(i,   j+h)) +
          C(pi<d>(i+1, j-h)) + 
          C(pi<d>(i,   j-h)) 
        ) / 4
      )

      template<int dim, class arr_2d_t>
      inline auto HOT_1_helper(
        const arr_2d_t &C,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
          (3 * C(pi<dim>(i+h, j)) * abs(C(pi<dim>(i+h, j))) - 2 * pow(C(pi<dim>(i+h, j)), 3) - C(pi<dim>(i+h, j))) / 3
      )

      template<int dim, class arr_2d_t>
      inline auto HOT_2_helper(
        const arrvec_t<arr_2d_t> &C,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
          (abs(C[dim](pi<dim>(i+h, j))) - 2 * pow(C[dim](pi<dim>(i+h, j)), 2)) / 4 * C_bar<dim>(C[dim-1], i, j)
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_1( // positive sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && opt_set(opts, pds)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+2, j)) - psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j))
            ,//-----------------------------------------------------------------------------------
            psi(pi<dim>(i+2, j)) + psi(pi<dim>(i+1, j)) + psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_1( // variable sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+2, j))) - abs(psi(pi<dim>(i+1, j))) - abs(psi(pi<dim>(i, j))) + abs(psi(pi<dim>(i-1, j)))
            ,//-------------------------------------------------------------------------------------------------------
            abs(psi(pi<dim>(i+2, j))) + abs(psi(pi<dim>(i+1, j))) + abs(psi(pi<dim>(i, j))) + abs(psi(pi<dim>(i-1, j)))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_1( // inf. gauge option
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opt_set(opts, iga)>::type* = 0
      ) return_macro(,
          ( psi(pi<dim>(i+2, j)) - psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j)) )
          / //--------------------------------------------------------------------------------------
          (1 + 1 + 1 + 1)
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_2( // positive sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && opt_set(opts, pds)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+1, j+1)) - psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1))            
            ,//-----------------------------------------------------------------------------------------
            psi(pi<dim>(i+1, j+1)) + psi(pi<dim>(i, j+1)) + psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_2( // variable sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+1, j+1))) - abs(psi(pi<dim>(i, j+1))) - abs(psi(pi<dim>(i+1, j-1))) + abs(psi(pi<dim>(i, j-1)))            
            ,//-------------------------------------------------------------------------------------------------------------
            abs(psi(pi<dim>(i+1, j+1))) + abs(psi(pi<dim>(i, j+1))) + abs(psi(pi<dim>(i+1, j-1))) + abs(psi(pi<dim>(i, j-1)))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_2( // inf. gauge option
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opt_set(opts, iga)>::type* = 0
      ) return_macro(,
          (psi(pi<dim>(i+1, j+1)) - psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1)))            
          / //-------------------------------------------------------------------------------------------
          (1 + 1 + 1 + 1)
      )

      // 3rd order term (see eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998 (with G=1))
      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT( // no HOT
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &C,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, toa)>::type* = 0
      ) -> decltype(0)
      {
        return 0; //no Higher Order Terms for second accuracy 
      }

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT( // higher order terms
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &C,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opt_set(opts, toa)>::type* = 0
      ) return_macro(,
         HOT_1<opts BOOST_PP_COMMA() dim>(psi, i, j) * HOT_1_helper<dim>(C[dim], i, j) 
         + 
         HOT_2<opts BOOST_PP_COMMA() dim>(psi, i, j) * HOT_2_helper<dim>(C, i, j)
      )

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto antidiff(
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        const arrvec_t<arr_2d_t> &C
      ) return_macro(,
        // second order terms
        abs(C[dim](pi<dim>(i+h, j))) 
        * (1 - abs(C[dim](pi<dim>(i+h, j)))) 
        * A<opts BOOST_PP_COMMA() dim>(psi, i, j) 
        - C[dim](pi<dim>(i+h, j)) 
        * C_bar<dim>(C[dim-1], i, j)
        * B<opts BOOST_PP_COMMA() dim>(psi, i, j)
        // third order terms
        + HOT<opts BOOST_PP_COMMA() dim>(psi, C, i, j)
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx 
