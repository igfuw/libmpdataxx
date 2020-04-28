/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_psi_2d.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      // divergent flow correction - no correction
      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto DFL(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::dfl)>::type* = 0
      )
      {
        return 0;
      }

      // divergent flow correction - general case
      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto DFL(
        const arr_2d_t &psi,    //to have the same arguments as in iga option
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::dfl) && !opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          - fconst<arr_2d_t>(0.25) * GC[dim](pi<dim>(i+h, j))
          /
          G_bar_x<opts, dim>(G, i, j)
          *
          (
            (
              GC[dim](pi<dim>((i+1)+h, j)) -
              GC[dim](pi<dim>(i-h    , j))
            )
            +
            (
              GC[dim-1](pi<dim>(i+1, j+h)) +
              GC[dim-1](pi<dim>(i,   j+h)) -
              GC[dim-1](pi<dim>(i+1, j-h)) -
              GC[dim-1](pi<dim>(i,   j-h))
            )
          )
        );
      }

      // divergent flow correction - infinite gauge version
      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      inline auto DFL(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::dfl) && opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          - fconst<arr_2d_t>(0.25) * GC[dim](pi<dim>(i+h, j))
          /
          G_bar_x<opts, dim>(G, i, j)
          *
          (
            (
              GC[dim](pi<dim>((i+1)+h, j)) -
              GC[dim](pi<dim>(i-h    , j))
            )
            +
            (
              GC[dim-1](pi<dim>(i+1, j+h)) +
              GC[dim-1](pi<dim>(i,   j+h)) -
              GC[dim-1](pi<dim>(i+1, j-h)) -
              GC[dim-1](pi<dim>(i,   j-h))
            )
          )
          * psi_bar_x<opts, dim>(psi, i, j)
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx
