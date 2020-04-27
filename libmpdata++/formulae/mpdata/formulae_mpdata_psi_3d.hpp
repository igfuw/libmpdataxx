/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the scalar field psi
// note that the function naming convention and the comments in this file correspond to dim = 0,
// if dim = 1 or 2 then you have to mentally perform appropiate substitutions
// for example for dim = 1 (x, y, z) -> (y, z, x) and (i+1/2, j, k) -> (i, j+1/2, k),
//             for dim = 2 (x, y, z) -> (z, x, y) and (i+1/2, j, k) -> (i, j, k+1/2)

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      // interpolation of psi to (i+1/2, j, k) - positive sign scalar / infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_x(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i+1, j, k)) +
            psi(pi<dim>(i  , j, k))
          ) / 2
        );
      }

      // interpolation of psi to (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_x(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i+1, j, k))) +
            abs(psi(pi<dim>(i  , j, k)))
          ) / 2
        );
      }

      // interpolation of psi to (i, j+1/2, k) - positive sign scalar / infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_y(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i, j  , k)) +
            psi(pi<dim>(i, j+1, k))
          ) / 2
        );
      }

      // interpolation of psi to (i, j+1/2, k) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_y(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i, j  , k))) +
            abs(psi(pi<dim>(i, j+1, k)))
          ) / 2
        );
      }

      // interpolation of psi to (i, j, k+1/2) - positive sign scalar / infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_z(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i, j, k  )) +
            psi(pi<dim>(i, j, k+1))
          ) / 2
        );
      }

      // interpolation of psi to (i, j, k+1/2) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_z(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i, j, k  ))) +
            abs(psi(pi<dim>(i, j, k+1)))
          ) / 2
        );
      }

      // interpolation of psi to (i+1/2, j+1/2, k) - positive sign scalar / infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_xy(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i+1, j  , k)) +
            psi(pi<dim>(i  , j  , k)) +
            psi(pi<dim>(i  , j+1, k)) +
            psi(pi<dim>(i+1, j+1, k))
          ) / 4
        );
      }

      // interpolation of psi to (i+1/2, j+1/2, k) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_xy(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i+1, j  , k))) +
            abs(psi(pi<dim>(i  , j  , k))) +
            abs(psi(pi<dim>(i  , j+1, k))) +
            abs(psi(pi<dim>(i+1, j+1, k)))
          ) / 4
        );
      }

      // interpolation of psi to (i+1/2, j, k+1/2) - positive sign scalar / infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_xz(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i+1, j  , k)) +
            psi(pi<dim>(i  , j  , k)) +
            psi(pi<dim>(i  , j, k+1)) +
            psi(pi<dim>(i+1, j, k+1))
          ) / 4
        );
      }

      // interpolation of psi to (i+1/2, j, k+1/2) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_xz(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i+1, j  , k))) +
            abs(psi(pi<dim>(i  , j  , k))) +
            abs(psi(pi<dim>(i  , j, k+1))) +
            abs(psi(pi<dim>(i+1, j, k+1)))
          ) / 4
        );
      }

      // interpolation of psi to (i+1/2, j+1/2, k+1/2) - positive sign scalar / infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_xyz(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            psi(pi<dim>(i  , j  , k  )) +
            psi(pi<dim>(i+1, j  , k  )) +
            psi(pi<dim>(i  , j+1, k  )) +
            psi(pi<dim>(i  , j  , k+1)) +
            psi(pi<dim>(i+1, j+1, k  )) +
            psi(pi<dim>(i+1, j  , k+1)) +
            psi(pi<dim>(i  , j+1, k+1)) +
            psi(pi<dim>(i+1, j+1, k+1))
          ) / 8
        );
      }

      // interpolation of psi to (i+1/2, j+1/2, k+1/2) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto psi_bar_xyz(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            abs(psi(pi<dim>(i  , j  , k  ))) +
            abs(psi(pi<dim>(i+1, j  , k  ))) +
            abs(psi(pi<dim>(i  , j+1, k  ))) +
            abs(psi(pi<dim>(i  , j  , k+1))) +
            abs(psi(pi<dim>(i+1, j+1, k  ))) +
            abs(psi(pi<dim>(i+1, j  , k+1))) +
            abs(psi(pi<dim>(i  , j+1, k+1))) +
            abs(psi(pi<dim>(i+1, j+1, k+1)))
          ) / 8
        );
      }

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndx_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
            psi(pi<d>(i+1, j, k)) - psi(pi<d>(i, j, k))
            ,
            psi(pi<d>(i+1, j, k)) + psi(pi<d>(i, j, k))
          )
        );
      }

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2, j, k) - variable-sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndx_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
            abs(psi(pi<d>(i+1, j, k))) - abs(psi(pi<d>(i, j, k)))
            ,
            abs(psi(pi<d>(i+1, j, k))) + abs(psi(pi<d>(i, j, k)))
          )
        );
      }

      // nondimensionalised x derivative of psi i.e.
      // dx/psi * dpsi/dx at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndx_psi(  // inf. gauge option
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        return return_helper<ix_t>(
          psi(pi<d>(i+1, j, k)) - psi(pi<d>(i, j, k))
        );
      }

      // nondimensionalised y derivative of psi i.e.
      // dy/psi * dpsi/dy at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndy_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          frac<opts, ix_t>(
              psi(pi<d>(i+1, j+1, k))
            + psi(pi<d>(i  , j+1, k))
            - psi(pi<d>(i+1, j-1, k))
            - psi(pi<d>(i  , j-1, k))
            ,
              psi(pi<d>(i+1, j+1, k))
            + psi(pi<d>(i  , j+1, k))
            + psi(pi<d>(i+1, j-1, k))
            + psi(pi<d>(i  , j-1, k))
          )
        );
      }

      // nondimensionalised y derivative of psi i.e.
      // dy/psi * dpsi/dy at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndy_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          frac<opts, ix_t>(
              abs(psi(pi<d>(i+1, j+1, k)))
            + abs(psi(pi<d>(i  , j+1, k)))
            - abs(psi(pi<d>(i+1, j-1, k)))
            - abs(psi(pi<d>(i  , j-1, k)))
            ,
              abs(psi(pi<d>(i+1, j+1, k)))
            + abs(psi(pi<d>(i  , j+1, k)))
            + abs(psi(pi<d>(i+1, j-1, k)))
            + abs(psi(pi<d>(i  , j-1, k)))
          )
        );
      }

      // nondimensionalised y derivative of psi i.e.
      // dy/psi * dpsi/dy at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndy_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        return return_helper<ix_t>(
          (
            psi(pi<d>(i+1, j+1, k))
          + psi(pi<d>(i  , j+1, k))
          - psi(pi<d>(i+1, j-1, k))
          - psi(pi<d>(i  , j-1, k))
          ) / 4
        );
      }

      // nondimensionalised z derivative of psi i.e.
      // dz/psi * dpsi/dz at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndz_psi( // positive sign signal
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          frac<opts, ix_t>(
              psi(pi<d>(i+1, j, k+1))
            + psi(pi<d>(i  , j, k+1))
            - psi(pi<d>(i+1, j, k-1))
            - psi(pi<d>(i  , j, k-1))
            ,
              psi(pi<d>(i+1, j, k+1))
            + psi(pi<d>(i  , j, k+1))
            + psi(pi<d>(i+1, j, k-1))
            + psi(pi<d>(i  , j, k-1))
          )
        );
      }

      // nondimensionalised z derivative of psi i.e.
      // dz/psi * dpsi/dz at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          frac<opts, ix_t>(
              abs(psi(pi<d>(i+1, j, k+1)))
            + abs(psi(pi<d>(i  , j, k+1)))
            - abs(psi(pi<d>(i+1, j, k-1)))
            - abs(psi(pi<d>(i  , j, k-1)))
            ,
              abs(psi(pi<d>(i+1, j, k+1)))
            + abs(psi(pi<d>(i  , j, k+1)))
            + abs(psi(pi<d>(i+1, j, k-1)))
            + abs(psi(pi<d>(i  , j, k-1)))
          )
        );
      }

      // nondimensionalised y derivative of psi i.e.
      // dz/psi * dpsi/dz at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        return return_helper<ix_t>(
          (
            psi(pi<d>(i+1, j, k+1))
          + psi(pi<d>(i  , j, k+1))
          - psi(pi<d>(i+1, j, k-1))
          - psi(pi<d>(i  , j, k-1))
          ) / 4
        );
      }

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxx_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
                psi(pi<dim>(i+2, j, k))
              - psi(pi<dim>(i+1, j, k))
              - psi(pi<dim>(i  , j, k))
              + psi(pi<dim>(i-1, j, k))
              ,
                psi(pi<dim>(i+2, j, k))
              + psi(pi<dim>(i+1, j, k))
              + psi(pi<dim>(i  , j, k))
              + psi(pi<dim>(i-1, j, k))
          )
        );
      }

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxx_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
                abs(psi(pi<dim>(i+2, j, k)))
              - abs(psi(pi<dim>(i+1, j, k)))
              - abs(psi(pi<dim>(i  , j, k)))
              + abs(psi(pi<dim>(i-1, j, k)))
              ,
                abs(psi(pi<dim>(i+2, j, k)))
              + abs(psi(pi<dim>(i+1, j, k)))
              + abs(psi(pi<dim>(i  , j, k)))
              + abs(psi(pi<dim>(i-1, j, k)))
          )
        );
      }

      // nondimensionalised xx derivative of psi i.e.
      // dx^2/psi * dpsi/dxx at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxx_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "iga & abs are mutually exclusive");
        return return_helper<ix_t>(
          2 * (
                psi(pi<dim>(i+2, j, k))
              - psi(pi<dim>(i+1, j, k))
              - psi(pi<dim>(i  , j, k))
              + psi(pi<dim>(i-1, j, k)) )
          / 4
        );
      }

      // nondimensionalised xy derivative of psi i.e.
      // dx*dy/psi * dpsi/dxdy at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxy_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
                psi(pi<dim>(i+1, j+1, k))
              - psi(pi<dim>(i  , j+1, k))
              - psi(pi<dim>(i+1, j-1, k))
              + psi(pi<dim>(i  , j-1, k))
              ,
                psi(pi<dim>(i+1, j+1, k))
              + psi(pi<dim>(i  , j+1, k))
              + psi(pi<dim>(i+1, j-1, k))
              + psi(pi<dim>(i  , j-1, k))
          )
        );
      }

      // nondimensionalised xy derivative of psi i.e.
      // dx*dy/psi * dpsi/dxdy at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxy_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
                abs(psi(pi<dim>(i+1, j+1, k)))
              - abs(psi(pi<dim>(i  , j+1, k)))
              - abs(psi(pi<dim>(i+1, j-1, k)))
              + abs(psi(pi<dim>(i  , j-1, k)))
              ,
                abs(psi(pi<dim>(i+1, j+1, k)))
              + abs(psi(pi<dim>(i  , j+1, k)))
              + abs(psi(pi<dim>(i+1, j-1, k)))
              + abs(psi(pi<dim>(i  , j-1, k)))
          )
        );
      }

      // nondimensionalised xy derivative of psi i.e.
      // dx*dy/psi * dpsi/dxdy at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxy_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        return return_helper<ix_t>(
          2 * (
                psi(pi<dim>(i+1, j+1, k))
              - psi(pi<dim>(i  , j+1, k))
              - psi(pi<dim>(i+1, j-1, k))
              + psi(pi<dim>(i  , j-1, k))
          ) / 4
        );
      }

      // nondimensionalised xz derivative of psi i.e.
      // dx*dz/psi * dpsi/dxdz at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
                psi(pi<dim>(i+1, j, k+1))
              - psi(pi<dim>(i  , j, k+1))
              - psi(pi<dim>(i+1, j, k-1))
              + psi(pi<dim>(i  , j, k-1))
              ,
                psi(pi<dim>(i+1, j, k+1))
              + psi(pi<dim>(i  , j, k+1))
              + psi(pi<dim>(i+1, j, k-1))
              + psi(pi<dim>(i  , j, k-1))
          )
        );
      }

      // nondimensionalised xz derivative of psi i.e.
      // dx*dz/psi * dpsi/dxdz at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            2 *
            frac<opts, ix_t>(
                abs(psi(pi<dim>(i+1, j, k+1)))
              - abs(psi(pi<dim>(i  , j, k+1)))
              - abs(psi(pi<dim>(i+1, j, k-1)))
              + abs(psi(pi<dim>(i  , j, k-1)))
              ,
                abs(psi(pi<dim>(i+1, j, k+1)))
              + abs(psi(pi<dim>(i  , j, k+1)))
              + abs(psi(pi<dim>(i+1, j, k-1)))
              + abs(psi(pi<dim>(i  , j, k-1)))
          )
        );
      }

      // nondimensionalised xz derivative of psi i.e.
      // dx*dz/psi * dpsi/dxdz at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndxz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        return return_helper<ix_t>(
          2 * (
                psi(pi<dim>(i+1, j, k+1))
              - psi(pi<dim>(i  , j, k+1))
              - psi(pi<dim>(i+1, j, k-1))
              + psi(pi<dim>(i  , j, k-1))
          ) / 4
        );
      }

      // nondimensionalised yz derivative of psi i.e.
      // dy*dz/psi * dpsi/dydz at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndyz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            frac<opts, ix_t>(
              psi(pi<dim>(i+1, j+1, k+1))
            + psi(pi<dim>(i+1, j-1, k-1))
            - psi(pi<dim>(i+1, j+1, k-1))
            - psi(pi<dim>(i+1, j-1, k+1))
            + psi(pi<dim>(i  , j+1, k+1))
            + psi(pi<dim>(i  , j-1, k-1))
            - psi(pi<dim>(i  , j+1, k-1))
            - psi(pi<dim>(i  , j-1, k+1))
            , //-------------------------
              psi(pi<dim>(i+1, j+1, k+1))
            + psi(pi<dim>(i+1, j-1, k-1))
            + psi(pi<dim>(i+1, j+1, k-1))
            + psi(pi<dim>(i+1, j-1, k+1))
            + psi(pi<dim>(i  , j+1, k+1))
            + psi(pi<dim>(i  , j-1, k-1))
            + psi(pi<dim>(i  , j+1, k-1))
            + psi(pi<dim>(i  , j-1, k+1))
            )
        );
      }

      // nondimensionalised yz derivative of psi i.e.
      // dy*dz/psi * dpsi/dydz at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndyz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            frac<opts, ix_t>(
              abs(psi(pi<dim>(i+1, j+1, k+1)))
            + abs(psi(pi<dim>(i+1, j-1, k-1)))
            - abs(psi(pi<dim>(i+1, j+1, k-1)))
            - abs(psi(pi<dim>(i+1, j-1, k+1)))
            + abs(psi(pi<dim>(i  , j+1, k+1)))
            + abs(psi(pi<dim>(i  , j-1, k-1)))
            - abs(psi(pi<dim>(i  , j+1, k-1)))
            - abs(psi(pi<dim>(i  , j-1, k+1)))
            , //------------------------------
              abs(psi(pi<dim>(i+1, j+1, k+1)))
            + abs(psi(pi<dim>(i+1, j-1, k-1)))
            + abs(psi(pi<dim>(i+1, j+1, k-1)))
            + abs(psi(pi<dim>(i+1, j-1, k+1)))
            + abs(psi(pi<dim>(i  , j+1, k+1)))
            + abs(psi(pi<dim>(i  , j-1, k-1)))
            + abs(psi(pi<dim>(i  , j+1, k-1)))
            + abs(psi(pi<dim>(i  , j-1, k+1)))
            )
        );
      }

      // nondimensionalised yz derivative of psi i.e.
      // dy*dz/psi * dpsi/dydz at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto ndyz_psi(
        const arr_3d_t &psi,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        return return_helper<ix_t>(
            (
              psi(pi<dim>(i+1, j+1, k+1))
            + psi(pi<dim>(i+1, j-1, k-1))
            - psi(pi<dim>(i+1, j+1, k-1))
            - psi(pi<dim>(i+1, j-1, k+1))
            + psi(pi<dim>(i  , j+1, k+1))
            + psi(pi<dim>(i  , j-1, k-1))
            - psi(pi<dim>(i  , j+1, k-1))
            - psi(pi<dim>(i  , j-1, k+1))
            )
            / //---------------------------
            (
              1 + 1 + 1 + 1 + 1 + 1 + 1 + 1
            )
        );
      }

      // nondimensionalised tx derivative of psi i.e.
      // dx*dt/psi * dpsi/dtx at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndtx_psi(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 *
          frac<opts, ix_t>(
              psi_np1(pi<d>(i+1, j, k))
            - psi_n(pi<d>(i+1, j, k))
            - psi_np1(pi<d>(i  , j, k))
            + psi_n(pi<d>(i  , j, k))
            ,
              psi_np1(pi<d>(i+1, j, k))
            + psi_n(pi<d>(i+1, j, k))
            + psi_np1(pi<d>(i  , j, k))
            + psi_n(pi<d>(i  , j, k))
          )
        );
      }

      // nondimensionalised tx derivative of psi i.e.
      // dx*dt/psi * dpsi/dtx at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndtx_psi(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 *
          frac<opts, ix_t>(
              abs(psi_np1(pi<d>(i+1, j, k)))
            - abs(psi_n(pi<d>(i+1, j, k)))
            - abs(psi_np1(pi<d>(i  , j, k)))
            + abs(psi_n(pi<d>(i  , j, k)))
            ,
              abs(psi_np1(pi<d>(i+1, j, k)))
            + abs(psi_n(pi<d>(i+1, j, k)))
            + abs(psi_np1(pi<d>(i  , j, k)))
            + abs(psi_n(pi<d>(i  , j, k)))
          )
        );
      }

      // nondimensionalised tx derivative of psi i.e.
      // dx*dt/psi * dpsi/dtx at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndtx_psi(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        return return_helper<ix_t>(
          psi_np1(pi<d>(i+1, j, k)) - psi_n(pi<d>(i+1, j, k)) - psi_np1(pi<d>(i, j, k)) + psi_n(pi<d>(i, j, k))
        );
      }

      // nondimensionalised t derivative of psi i.e.
      // dt/psi * dpsi/dt at (i+1/2, j, k) - positive sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndt_psi(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
              psi_np1(pi<d>(i+1, j, k))
            - psi_n(pi<d>(i+1, j, k))
            + psi_np1(pi<d>(i  , j, k))
            - psi_n(pi<d>(i  , j, k))
            ,
              psi_np1(pi<d>(i+1, j, k))
            + psi_n(pi<d>(i+1, j, k))
            + psi_np1(pi<d>(i  , j, k))
            + psi_n(pi<d>(i  , j, k))
          )
        );
      }

      // nondimensionalised t derivative of psi i.e.
      // dt/psi * dpsi/dt at (i+1/2, j, k) - variable sign scalar version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndt_psi(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
              abs(psi_np1(pi<d>(i+1, j, k)))
            - abs(psi_n(pi<d>(i+1, j, k)))
            + abs(psi_np1(pi<d>(i  , j, k)))
            - abs(psi_n(pi<d>(i  , j, k)))
            ,
              abs(psi_np1(pi<d>(i+1, j, k)))
            + abs(psi_n(pi<d>(i+1, j, k)))
            + abs(psi_np1(pi<d>(i  , j, k)))
            + abs(psi_n(pi<d>(i  , j, k)))
          )
        );
      }

      // nondimensionalised tx derivative of psi i.e.
      // dt/psi * dpsi/dt at (i+1/2, j, k) - infinite gauge version
      template<opts_t opts, int d, class arr_3d_t, class ix_t>
      inline auto ndt_psi(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        return return_helper<ix_t>(
          fconst<arr_3d_t>(0.5) * (psi_np1(pi<d>(i+1, j, k)) - psi_n(pi<d>(i+1, j, k)) + psi_np1(pi<d>(i, j, k)) - psi_n(pi<d>(i, j, k)))
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
