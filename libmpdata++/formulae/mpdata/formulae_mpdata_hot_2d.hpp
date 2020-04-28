/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_psi_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_gc_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_2d.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto ndxx_psi_coeff(
        const arr_2d_t &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          (
            3 * GC(pi<dim>(i+h, j)) * abs(GC(pi<dim>(i+h, j))) / G_bar_x<opts, dim>(G, i, j)
            - 2 * pow3(GC(pi<dim>(i+h, j))) / pow2(G_bar_x<opts, dim>(G, i, j))
            - GC(pi<dim>(i+h, j))
          ) / 6
        );
      }

      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto ndxy_psi_coeff(
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          (abs(GC[dim](pi<dim>(i+h, j))) - 2 * pow2(GC[dim](pi<dim>(i+h, j))) / G_bar_x<opts, dim>(G, i, j))
           * GC1_bar_xy<dim>(GC[dim+1], i, j) / (2 * G_bar_x<opts, dim>(G, i, j))
        );
      }

      // third order terms
      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto TOT(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::tot)>::type* = 0
      )
      {
        return return_helper<ix_t>(
           ndxx_psi<opts, dim>(psi, i, j) * ndxx_psi_coeff<opts, dim>(GC[dim], G, i, j)
           +
           ndxy_psi<opts, dim>(psi, i, j) * ndxy_psi_coeff<opts, dim>(GC, G, i, j)
        );
      }

      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto TOT(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::tot)>::type* = 0
      )
      {
        return 0;
      }

      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto ndxxx_psi_coeff(
        const arr_2d_t &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          (
            6 * abs(GC(pi<dim>(i+h, j))) * pow2(GC(pi<dim>(i+h, j))) / pow2(G_bar_x<opts, dim>(G, i, j))
            - 3 * pow2(GC(pi<dim>(i+h, j))) / G_bar_x<opts, dim>(G, i, j)
            - 3 * pow4(GC(pi<dim>(i+h, j))) / pow3(G_bar_x<opts, dim>(G, i, j))
          ) / 2
        );
      }

      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto ndxxy_psi_coeff(
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          3 *
          (abs(GC[dim](pi<dim>(i+h, j))) - pow2(GC[dim](pi<dim>(i+h, j))) / G_bar_x<opts, dim>(G, i, j))
           * GC[dim](pi<dim>(i+h, j)) * GC1_bar_xy<dim>(GC[dim-1], i, j)
           / pow2(G_bar_x<opts, dim>(G, i, j))
        );
      }

      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto ndxyy_psi_coeff(
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          3 *
          (
              3 * abs(GC1_bar_xy<dim>(GC[dim-1], i, j)) * pow2(GC[dim](pi<dim>(i+h, j))) / G_bar_x<opts, dim>(G, i, j)
            + 3 * pow2(GC1_bar_xy<dim>(GC[dim-1], i, j)) * abs(GC[dim](pi<dim>(i+h, j))) / G_bar_x<opts, dim>(G, i, j)
            - 2 * abs(GC1_bar_xy<dim>(GC[dim-1], i, j)) * abs(GC[dim](pi<dim>(i+h, j)))
            - 9 * pow2(GC1_bar_xy<dim>(GC[dim-1], i, j)) * pow2(GC[dim](pi<dim>(i+h, j)))
                / pow2(G_bar_x<opts, dim>(G, i, j))
          ) / G_bar_x<opts, dim>(G, i, j) / 4
        );
      }

      // fourth order terms
      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto FOT(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::fot)>::type* = 0
      )
      {
        static_assert(opts::isset(opts, opts::tot) || opts::isset(opts, opts::div_3rd),
                      "adding fourth-order terms makes sense only when third-order terms are present (tot or div_3rd option)");
        static_assert(opts::isset(opts, opts::iga), "fot option only available with iga");
        return return_helper<ix_t>(
         ndxxx_psi<opts, dim>(psi, i, j) * ndxxx_psi_coeff<opts, dim>(GC[dim], G, i, j)
         +
         ndxxy_psi<opts, dim>(psi, i, j) * ndxxy_psi_coeff<opts, dim>(GC, G, i, j)
         +
         ndxyy_psi<opts, dim>(psi, i, j) * ndxyy_psi_coeff<opts, dim>(GC, G, i, j)
        );
      }

      template<opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto FOT(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::fot)>::type* = 0
      )
      {
        return 0;
      }
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
