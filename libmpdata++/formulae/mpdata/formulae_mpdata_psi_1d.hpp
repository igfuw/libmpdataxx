/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the scalar field psi

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // interpolation of psi to (i+1/2) - positive sign scalar / infinite gauge version
      template<opts_t opts, class arr_1d_t>
      inline auto psi_bar_x( 
        const arr_1d_t &psi,
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        (
          psi(i+1) + psi(i)
        ) / 2
      )
     
      // interpolation of psi to (i+1/2) - variable sign scalar version
      template<opts_t opts, class arr_1d_t>
      inline auto psi_bar_x( 
        const arr_1d_t &psi,
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        (
          abs(psi(i+1)) + abs(psi(i))
        ) / 2
      )

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2) - positive sign scalar version  
      template<opts_t opts, class arr_1d_t>
      inline auto ndx_psi(
        const arr_1d_t &psi, 
        const rng_t &i, 
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        2 *
	frac<opts>(
	  psi(i+1) - psi(i)
	  ,// -------------
	  psi(i+1) + psi(i)
	)
      )

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2) - variable-sign scalar version 
      template<opts_t opts, class arr_1d_t>
      inline auto ndx_psi(
        const arr_1d_t &psi, 
        const rng_t &i, 
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        2 *
	frac<opts>(
	  abs(psi(i+1)) - abs(psi(i))
	  ,// -----------------------
	  abs(psi(i+1)) + abs(psi(i))
	) 
      ) 

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2) - infinite gauge version 
      template<opts_t opts, class arr_1d_t>
      inline auto ndx_psi(
        const arr_1d_t &psi, 
        const rng_t &i, 
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        ,
        2 *
        (
	  psi(i+1) - psi(i)
        ) / ( //-----------
	  1 + 1
        )
      )
      
      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2) - positive sign scalar version
      template<opts_t opts, class arr_1d_t>
      inline auto ndxx_psi(
        const arr_1d_t &psi,
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          2 *
          frac<opts>(
            psi(i+2) - psi(i+1) - psi(i) + psi(i-1)
            ,//------------------------------------
            psi(i+2) + psi(i+1) + psi(i) + psi(i-1)
        )
      )

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2) - variable sign scalar version
      template<opts_t opts, class arr_1d_t>
      inline auto ndxx_psi(
        const arr_1d_t &psi,
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          2 *
          frac<opts>(
            abs(psi(i+2)) - abs(psi(i+1)) - abs(psi(i)) + abs(psi(i-1))
            ,//--------------------------------------------------------
            abs(psi(i+2)) + abs(psi(i+1)) + abs(psi(i)) + abs(psi(i-1))
        )
      )

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2) - variable sign scalar version
      template<opts_t opts, class arr_1d_t>
      inline auto ndxx_psi(
        const arr_1d_t &psi,
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs are mutually exclusive");
        ,
        2 *
	( psi(i+2) - psi(i+1) - psi(i) + psi(i-1) )
	/ //---------------------------------------
	(1 + 1 + 1 + 1)
      )
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx 
