/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the scalar field psi
// note that the function naming convention and the comments in this file correspond to dim = 0,
// if dim = 1 then you have to mentally perform appropiate substitutions
// for example (x, y) -> (y, x) and (i+1/2, j) -> (i, j+1/2)

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      // interpolation of psi to (i+1/2, j) - positive sign scalar / infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto psi_bar_x(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i+1, j)) + psi(pi<dim>(i, j))
          ) / 2
        );
      }

      // interpolation of psi to (i+1/2, j) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto psi_bar_x(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i+1, j))) + abs(psi(pi<dim>(i, j)))
          ) / 2
        );
      }

      // interpolation of psi to (i, j+1/2) - positive sign scalar / infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto psi_bar_y(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i, j)) + psi(pi<dim>(i, j+1))
          ) / 2
        );
      }

      // interpolation of psi to (i, j+1/2) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto psi_bar_y(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i, j))) + abs(psi(pi<dim>(i, j+1)))
          ) / 2
        );
      }

      // interpolation of psi to (i+1/2, j+1/2) - positive sign scalar / infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto psi_bar_xy(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i+1, j  )) +
            psi(pi<dim>(i  , j  )) +
            psi(pi<dim>(i  , j+1)) +
            psi(pi<dim>(i+1, j+1))
          ) / 4
        );
      }

      // interpolation of psi to (i+1/2, j+1/2) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto psi_bar_xy(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i+1, j  ))) +
            abs(psi(pi<dim>(i  , j  ))) +
            abs(psi(pi<dim>(i  , j+1))) +
            abs(psi(pi<dim>(i+1, j+1)))
          ) / 4
        );
      }

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2, j) - positive sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndx_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
            psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j))
            ,// ---------------------------------
            psi(pi<dim>(i+1, j)) + psi(pi<dim>(i, j))
          )
        );
      }

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2, j) - variable-sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndx_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
            abs(psi(pi<dim>(i+1, j))) - abs(psi(pi<dim>(i, j)))
            ,// -------------------------------------------
            abs(psi(pi<dim>(i+1, j))) + abs(psi(pi<dim>(i, j)))
          )
        );
      }

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndx_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        return return_helper<ix_t>(
          2 *
          (
            psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j))
          ) / ( //-------------------------------
            1 + 1
          )
        );
      }

      // nondimensionalised y derivative of psi i.e.
      // dy/psi * dpsi/dy at (i+1/2, j) - positive sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndy_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          frac<opts, ix_t>(
            psi(pi<dim>(i+1, j+1)) + psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) - psi(pi<dim>(i, j-1))
            ,// ---------------------------------------------------------------------------------
            psi(pi<dim>(i+1, j+1)) + psi(pi<dim>(i, j+1)) + psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1))
          )
        );
      }

      // nondimensionalised y derivative of psi i.e.
      // dy/psi * dpsi/dy at (i+1/2, j) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndy_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          frac<opts, ix_t>(
            abs(psi(pi<dim>(i+1, j+1))) + abs(psi(pi<dim>(i, j+1))) - abs(psi(pi<dim>(i+1, j-1))) - abs(psi(pi<dim>(i, j-1)))
            ,// -----------------------------------------------------------------------------------------------------
            abs(psi(pi<dim>(i+1, j+1))) + abs(psi(pi<dim>(i, j+1))) + abs(psi(pi<dim>(i+1, j-1))) + abs(psi(pi<dim>(i, j-1)))
          )
        );
      }

      // nondimensionalised y derivative of psi i.e.
      // dy/psi * dpsi/dy at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndy_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i+1, j+1)) + psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) - psi(pi<dim>(i, j-1))
          ) / (  // -----------------------------------------------------------------------------
            1 + 1 + 1 + 1
          )
        );
      }

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2, j) - positive sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxx_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
              psi(pi<dim>(i+2, j)) - psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j))
              ,//------------------------------------------------------------------------------------
              psi(pi<dim>(i+2, j)) + psi(pi<dim>(i+1, j)) + psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j))
          )
        );
      }

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2, j) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxx_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
              abs(psi(pi<dim>(i+2, j))) - abs(psi(pi<dim>(i+1, j))) - abs(psi(pi<dim>(i, j))) + abs(psi(pi<dim>(i-1, j)))
              ,//-------------------------------------------------------------------------------------------------------
              abs(psi(pi<dim>(i+2, j))) + abs(psi(pi<dim>(i+1, j))) + abs(psi(pi<dim>(i, j))) + abs(psi(pi<dim>(i-1, j)))
          )
        );
      }

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxx_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "iga & abs are mutually exclusive");
        return return_helper<ix_t>(
          2 *
          ( psi(pi<dim>(i+2, j)) - psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j)) )
          / //--------------------------------------------------------------------------------------
          (1 + 1 + 1 + 1)
        );
      }

      // nondimensionalised xy derivative of psi i.e.
      // dx*dy/psi * dpsi/dxdy at (i+1/2, j) - positive sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxy_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
              psi(pi<dim>(i+1, j+1)) - psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1))
              ,//-----------------------------------------------------------------------------------------
              psi(pi<dim>(i+1, j+1)) + psi(pi<dim>(i, j+1)) + psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1))
          )
        );
      }

      // nondimensionalised xy derivative of psi i.e.
      // dx*dy/psi * dpsi/dxdy at (i+1/2, j) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxy_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
              abs(psi(pi<dim>(i+1, j+1))) - abs(psi(pi<dim>(i, j+1))) - abs(psi(pi<dim>(i+1, j-1))) + abs(psi(pi<dim>(i, j-1)))
              ,//-------------------------------------------------------------------------------------------------------------
              abs(psi(pi<dim>(i+1, j+1))) + abs(psi(pi<dim>(i, j+1))) + abs(psi(pi<dim>(i+1, j-1))) + abs(psi(pi<dim>(i, j-1)))
          )
        );
      }

      // nondimensionalised xy derivative of psi i.e.
      // dx*dy/psi * dpsi/dxdy at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxy_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        return return_helper<ix_t>(
          2 *
          (psi(pi<dim>(i+1, j+1)) - psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1)))
          / //-------------------------------------------------------------------------------------------
          (1 + 1 + 1 + 1)
        );
      }

      // nondimensionalised xxx derivative of psi i.e.
      // dx^3/psi * dpsi/dxxx at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxxx_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ( psi(pi<dim>(i+2, j)) - 3 * psi(pi<dim>(i+1, j)) + 3 * psi(pi<dim>(i, j)) - psi(pi<dim>(i-1, j)) )
          / //-----------------------------------------------------------------------------------------------
          (1 + 1 + 1 + 1)
        );
      }

      // nondimensionalised xxy derivative of psi i.e.
      // dx^2*dy/psi * dpsi/dxxy at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxxy_psi(
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (  psi(pi<dim>(i+2, j+1)) - psi(pi<dim>(i+1, j+1)) - psi(pi<dim>(i, j+1)) + psi(pi<dim>(i-1, j+1))
            -psi(pi<dim>(i+2, j-1)) + psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1)) - psi(pi<dim>(i-1, j-1)) )
          / //------------------------------------------------------------------------------------------------
          (1 + 1 + 1 + 1 + 1 + 1 + 1 + 1)
        );
      }

      // nondimensionalised xyy derivative of psi i.e.
      // dx*dy^2/psi * dpsi/dxyy at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto ndxyy_psi( // inf. gauge option
        const arr_2d_t &psi,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ( psi(pi<dim>(i+1, j+1)) - 2 * psi(pi<dim>(i+1, j)) + psi(pi<dim>(i+1, j-1))
           -psi(pi<dim>(i  , j+1)) + 2 * psi(pi<dim>(i  , j)) - psi(pi<dim>(i  , j-1)) )
          / //--------------------------------------------------------------------------
          (1 + 1 + 1 + 1 + 1 + 1)
        );
      }

      // nondimensionalised tx derivative of psi i.e.
      // dx*dt/psi * dpsi/dtx at (i+1/2, j) - positive sign scalar version
      template<opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto ndtx_psi(
        const arr_2d_t &psi_np1,
        const arr_2d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 *
          frac<opts, ix_t>(
              psi_np1(pi<d>(i+1, j))
            - psi_n(pi<d>(i+1, j))
            - psi_np1(pi<d>(i  , j))
            + psi_n(pi<d>(i  , j))
            ,
              psi_np1(pi<d>(i+1, j))
            + psi_n(pi<d>(i+1, j))
            + psi_np1(pi<d>(i  , j))
            + psi_n(pi<d>(i  , j))
          )
        );
      }

      // nondimensionalised tx derivative of psi i.e.
      // dx*dt/psi * dpsi/dtx at (i+1/2, j) - variable sign scalar version
      template<opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto ndtx_psi(
        const arr_2d_t &psi_np1,
        const arr_2d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 *
          frac<opts, ix_t>(
              abs(psi_np1(pi<d>(i+1, j)))
            - abs(psi_n(pi<d>(i+1, j)))
            - abs(psi_np1(pi<d>(i  , j)))
            + abs(psi_n(pi<d>(i  , j)))
            ,
              abs(psi_np1(pi<d>(i+1, j)))
            + abs(psi_n(pi<d>(i+1, j)))
            + abs(psi_np1(pi<d>(i  , j)))
            + abs(psi_n(pi<d>(i  , j)))
          )
        );
      }

      // nondimensionalised tx derivative of psi i.e.
      // dx*dt/psi * dpsi/dtx at (i+1/2, j) - infinite gauge version
      template<opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto ndtx_psi(
        const arr_2d_t &psi_np1,
        const arr_2d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            psi_np1(pi<d>(i+1, j))
          - psi_n(pi<d>(i+1, j))
          - psi_np1(pi<d>(i  , j))
          + psi_n(pi<d>(i  , j))
        );
      }

      // nondimensionalised t derivative of psi i.e.
      // dt/psi * dpsi/dt at (i+1/2, j) - positive sign scalar version
      template<opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto ndt_psi(
        const arr_2d_t &psi_np1,
        const arr_2d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
              psi_np1(pi<d>(i+1, j))
            - psi_n(pi<d>(i+1, j))
            + psi_np1(pi<d>(i  , j))
            - psi_n(pi<d>(i  , j))
            ,
              psi_np1(pi<d>(i+1, j))
            + psi_n(pi<d>(i+1, j))
            + psi_np1(pi<d>(i  , j))
            + psi_n(pi<d>(i  , j))
          )
        );
      }

      // nondimensionalised t derivative of psi i.e.
      // dt/psi * dpsi/dt at (i+1/2, j) - variable sign scalar version
      template<opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto ndt_psi(
        const arr_2d_t &psi_np1,
        const arr_2d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
              abs(psi_np1(pi<d>(i+1, j)))
            - abs(psi_n(pi<d>(i+1, j)))
            + abs(psi_np1(pi<d>(i  , j)))
            - abs(psi_n(pi<d>(i  , j)))
            ,
              abs(psi_np1(pi<d>(i+1, j)))
            + abs(psi_n(pi<d>(i+1, j)))
            + abs(psi_np1(pi<d>(i  , j)))
            + abs(psi_n(pi<d>(i  , j)))
          )
        );
      }

      // nondimensionalised t derivative of psi i.e.
      // dt/psi * dpsi/dt at (i+1/2, j) - infinite gauge version
      template<opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto ndt_psi(
        const arr_2d_t &psi_np1,
        const arr_2d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi_np1(pi<d>(i+1, j))
          - psi_n(pi<d>(i+1, j))
          + psi_np1(pi<d>(i  , j))
          - psi_n(pi<d>(i  , j))
          ) / 2
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
