/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the advector field GC
// note that the function naming convention and the comments in this file correspond to dim = 0,
// if dim = 1 then you have to mentally perform appropiate substitutions
// for example (GC[0], GC[1]) -> (GC[1], GC[0]), (x, y) -> (y, x) and (i+1/2, j) -> (i, j+1/2)

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_2d.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // interpolation of GC[0] to (i, j)
      template <int dim, class arr_2d_t, class ix_t>
      inline auto GC0_bar_x( 
        const arr_2d_t &GC,
        const ix_t &i,
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+h, j)) + GC(pi<dim>(i-h, j))
          ) / 2
        );
      }

      // interpolation of GC[0] to (i, j+1/2)
      template <int dim, class arr_2d_t, class ix_t>
      inline auto GC0_bar_xy(
        const arr_2d_t &GC, 
        const ix_t &i, 
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+h, j  )) + 
            GC(pi<dim>(i-h, j  )) +
            GC(pi<dim>(i+h, j+1)) + 
            GC(pi<dim>(i-h, j+1)) 
          ) / 4
        );
      }
      
      
      // interpolation of GC[1] to (i+1/2, j+1/2)
      // caution proper call looks like GC1_bar_x<dim>(GC[dim+1], i, j) - note dim vs dim+1 
      template <int dim, class arr_2d_t, class ix_t>
      inline auto GC1_bar_x( 
        const arr_2d_t &GC,
        const ix_t &i,
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j+h)) + GC(pi<dim>(i, j+h))
          ) / 2
        );
      }
 
      // interpolation of GC[1] to (i+1/2, j)
      // caution proper call looks like GC1_bar_xy<dim>(GC[dim+1], i, j) - note dim vs dim+1 
      template <int dim, class arr_2d_t, class ix_t>
      inline auto GC1_bar_xy(
        const arr_2d_t &GC, 
        const ix_t &i, 
        const ix_t &j
      )
      { 
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j+h)) + 
            GC(pi<dim>(i,   j+h)) +
            GC(pi<dim>(i+1, j-h)) + 
            GC(pi<dim>(i,   j-h)) 
          ) / 4
        );
      }

      // nondimensionalised x derivative of GC[0] i.e.
      // dx * dGC[0]/dx at (i+1/2, j)
      template <int dim, class arr_2d_t, class ix_t>
      inline auto ndx_GC0(
        const arr_2d_t &GC, 
        const ix_t &i,
        const ix_t &j
      )
      { 
        return return_helper<ix_t>(
          (GC(pi<dim>(i+h+1, j)) - GC(pi<dim>(i+h-1, j))) / 2
        );
      }
      
      // nondimensionalised xx derivative of GC[0] i.e.
      // dx^2 * dGC[0]/dxx at (i+1/2, j) - general case
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxx_GC0(
        const arr_2d_t &psi, 
        const arr_2d_t &GC, 
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          GC(pi<dim>(i+h+1, j)) + GC(pi<dim>(i+h-1, j)) - 2 * GC(pi<dim>(i+h, j))
        );
      }
      
      // nondimensionalised xx derivative of GC[0] i.e.
      // dx^2 * dGC[0]/dxx at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxx_GC0(
        const arr_2d_t &psi, 
        const arr_2d_t &GC, 
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (GC(pi<dim>(i+h+1, j)) + GC(pi<dim>(i+h-1, j)) - 2 * GC(pi<dim>(i+h, j))) * 
           psi_bar_x<opts, dim>(psi, i, j)
        );
      }
      
      // nondimensionalised tt derivative of GC[0] i.e.
      // dt^2 * dGC[0]/dtt at (i+1/2, j) - general case
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndtt_GC0(
        const arr_2d_t &psi, 
        const arr_2d_t &ndtt_GC, 
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndtt_GC(pi<dim>(i+h, j)) + 0
        );
      }
      
      // nondimensionalised tt derivative of GC[0] i.e.
      // dt^2 * dGC[0]/dtt at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndtt_GC0(
        const arr_2d_t &psi, 
        const arr_2d_t &ndtt_GC, 
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndtt_GC(pi<dim>(i+h, j)) * psi_bar_x<opts, dim>(psi, i, j)
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
