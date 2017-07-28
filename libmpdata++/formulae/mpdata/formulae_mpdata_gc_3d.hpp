/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the advector field GC
// note that the function naming convention and the comments in this file correspond to dim = 0,
// if dim = 1 or 2 then you have to mentally perform appropiate substitutions
// for example for dim = 1 (GC[0], GC[1], GC[2]) -> (GC[1], GC[2], GC[0]), 
//                         (x, y, z) -> (y, z, x) and (i+1/2, j, k) -> (i, j+1/2, k),
//             for dim = 2 (GC[0], GC[1], GC[2]) -> (GC[2], GC[0], GC[1]), 
//                         (x, y, z) -> (z, x, y) and (i+1/2, j, k) -> (i, j, k+1/2)

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_2d.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // interpolation of GC[0] to (i, j, k)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC0_bar_x( 
        const arr_3d_t &GC,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+h, j, k)) + 
            GC(pi<dim>(i-h, j, k))
          ) / 2
        );
      }
      
      // interpolation of GC[0] to (i, j+1/2, k)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC0_bar_xy(
        const arr_3d_t &GC, 
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+h, j  , k)) + 
            GC(pi<dim>(i-h, j  , k)) +
            GC(pi<dim>(i+h, j+1, k)) + 
            GC(pi<dim>(i-h, j+1, k)) 
          ) / 4
        );
      }

      // interpolation of GC[0] to (i, j, k+1/2)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC0_bar_xz(
        const arr_3d_t &GC, 
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+h, j, k  )) + 
            GC(pi<dim>(i-h, j, k  )) +
            GC(pi<dim>(i+h, j, k+1)) + 
            GC(pi<dim>(i-h, j, k+1)) 
          ) / 4
        );
      }

      // interpolation of GC[1] to (i+1/2, j+1/2, k)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC1_bar_x( 
        const arr_3d_t &GC,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j+h, k)) + 
            GC(pi<dim>(i  , j+h, k))
          ) / 2
        );
      }

      // interpolation of GC[1] to (i+1/2, j, k)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC1_bar_xy(
        const arr_3d_t &GC, 
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j+h, k)) + 
            GC(pi<dim>(i  , j+h, k)) +
            GC(pi<dim>(i+1, j-h, k)) + 
            GC(pi<dim>(i  , j-h, k)) 
          ) / 4
        );
      }
      
      // interpolation of GC[1] to (i+1/2, j+1/2, k+1/2)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC1_bar_xz(
        const arr_3d_t &GC, 
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j+h, k  )) + 
            GC(pi<dim>(i+1, j+h, k+1)) +
            GC(pi<dim>(i  , j+h, k  )) + 
            GC(pi<dim>(i  , j+h, k+1)) 
          ) / 4
        );
      }
      
      // interpolation of GC[2] to (i+1/2, j, k)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC2_bar_x(
        const arr_3d_t &GC, 
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j, k+h)) + 
            GC(pi<dim>(i  , j, k+h))
          ) / 2
        );
      }

      // interpolation of GC[2] to (i+1/2, j, k)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC2_bar_xz(
        const arr_3d_t &GC, 
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j, k+h)) + 
            GC(pi<dim>(i  , j, k+h)) +
            GC(pi<dim>(i+1, j, k-h)) + 
            GC(pi<dim>(i  , j, k-h)) 
          ) / 4
        );
      }
      
      // interpolation of GC[2] to (i+1/2, j+1/2, k+1/2)
      template<int dim, class arr_3d_t, class ix_t>
      inline auto GC2_bar_xy(
        const arr_3d_t &GC, 
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC(pi<dim>(i+1, j+1, k+h)) + 
            GC(pi<dim>(i  , j+1, k+h)) +
            GC(pi<dim>(i+1, j  , k+h)) + 
            GC(pi<dim>(i  , j  , k+h))
          ) / 4
        );
      }

      // nondimensionalised x derivative of GC[0] i.e.
      // dx * dGC[0]/dx at (i+1/2, j, k)
      template <int dim, class arr_3d_t, class ix_t>
      inline auto ndx_GC0(
        const arr_3d_t &GC, 
        const ix_t &i,
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (GC(pi<dim>(i+h+1, j, k)) - GC(pi<dim>(i+h-1, j, k))) / 2
        );
      }
      
      // nondimensionalised xx derivative of GC[0] i.e.
      // dx^2 * dGC[0]/dxx at (i+1/2, j, k) - general case
      template <opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxx_GC0(
        const arr_3d_t &psi, 
        const arr_3d_t &GC, 
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          GC(pi<dim>(i+h+1, j, k)) + GC(pi<dim>(i+h-1, j, k)) - 2 * GC(pi<dim>(i+h, j, k))
        );
      }
      
      // nondimensionalised xx derivative of GC[0] i.e.
      // dx^2 * dGC[0]/dxx at (i+1/2, j, k) - infinite gauge version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxx_GC0(
        const arr_3d_t &psi, 
        const arr_3d_t &GC, 
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (GC(pi<dim>(i+h+1, j, k)) + GC(pi<dim>(i+h-1, j, k)) - 2 * GC(pi<dim>(i+h, j, k))) * 
           psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
        );
      }
      
      // nondimensionalised tt derivative of GC[0] i.e.
      // dt^2 * dGC[0]/dtt at (i+1/2, j, k) - general case
      template <opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndtt_GC0(
        const arr_3d_t &psi, 
        const arr_3d_t &ndtt_GC,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndtt_GC(pi<dim>(i+h, j, k)) + 0
        );
      }
      
      // nondimensionalised tt derivative of GC[0] i.e.
      // dt^2 * dGC[0]/dtt at (i+1/2, j, k) - infinite gauge version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndtt_GC0(
        const arr_3d_t &psi, 
        const arr_3d_t &ndtt_GC,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndtt_GC(pi<dim>(i+h, j, k)) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
