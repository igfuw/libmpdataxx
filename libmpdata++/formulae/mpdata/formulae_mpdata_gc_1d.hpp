/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the advector field GC

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_2d.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // interpolation of GC[0] to (i)
      template<class arr_1d_t, class ix_t>
      inline auto GC0_bar_x( 
        const arr_1d_t &GC,
        const ix_t &i
      )
      {
        return return_helper<ix_t>(
          (
            GC(i+h) + GC(i-h)
          ) / 2
        );
      }

      // nondimensionalised x derivative of GC[0] i.e.
      // dx * dGC[0]/dx at (i+1/2)
      template <class arr_1d_t, class ix_t>
      inline auto ndx_GC0(
        const arr_1d_t &GC, 
        const ix_t &i
      )
      {
        return return_helper<ix_t>(
          (GC(i+h+1) - GC(i+h-1)) / 2
        );
      }
      
      // nondimensionalised xx derivative of GC[0] i.e.
      // dx^2 * dGC[0]/dxx at (i+1/2) - general case
      template <opts_t opts, class arr_1d_t, class ix_t>
      inline auto ndxx_GC0(
        const arr_1d_t &psi, 
        const arr_1d_t &GC, 
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          GC(i+h+1) + GC(i+h-1) - 2 * GC(i+h)
        );
      }
      
      // nondimensionalised xx derivative of GC[0] i.e.
      // dx^2 * dGC[0]/dxx at (i+1/2) - infinite gauge version
      template <opts_t opts, class arr_1d_t, class ix_t>
      inline auto ndxx_GC0(
        const arr_1d_t &psi, 
        const arr_1d_t &GC, 
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (GC(i+h+1) + GC(i+h-1) - 2 * GC(i+h)) * 
           psi_bar_x<opts>(psi, i)
        );
      }
      
      // nondimensionalised tt derivative of GC[0] i.e.
      // dt^2 * dGC[0]/dtt at (i+1/2) - general case
      template <opts_t opts, class arr_1d_t, class ix_t>
      inline auto ndtt_GC0(
        const arr_1d_t &psi, 
        const arr_1d_t &ndtt_GC, 
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndtt_GC(i+h) + 0
        );
      }
      
      // nondimensionalised tt derivative of GC[0] i.e.
      // dGC[0]/dtt at (i+1/2) - infinite-gauge version
      template <opts_t opts, class arr_1d_t, class ix_t>
      inline auto ndtt_GC0(
        const arr_1d_t &psi, 
        const arr_1d_t &ndtt_GC, 
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndtt_GC(i+h) * psi_bar_x<opts>(psi, i)
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
