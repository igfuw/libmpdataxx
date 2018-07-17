/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_dfl_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_hot_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fdiv_2d.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // first come helpers for divergence form of antidiffusive velocity
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto div_2nd(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G, 
        const ix_t &i, 
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          // second order terms
          abs(GC[dim](pi<dim>(i+h, j))) / 2
          * ndx_psi<opts, dim>(psi, i, j) 
          - 
          GC[dim](pi<dim>(i+h, j)) / 2
          * nfdiv<opts, dim>(psi, GC, G, i, j)
        );
      }
      
      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_upwind(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G, 
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          abs(div_2nd<opts, dim>(psi, GC, G, i, j)) / 2
          * ndx_psi<opts, dim>(psi, i, j) 
        );
      }

      template <opts_t opts, int dim, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_upwind(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G, 
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return 0;
      }
      
      template <opts_t opts, int dim, solvers::tmprl_extrp_t tmprl_extrp, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_temporal(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &ndtt_GC,
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<tmprl_extrp == solvers::noextrp>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndtt_GC0<opts, dim>(psi, ndtt_GC[dim], i, j)
        );
      }
      
      template <opts_t opts, int dim, solvers::tmprl_extrp_t tmprl_extrp, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_temporal(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &ndtt_GC,
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<tmprl_extrp == solvers::linear2>::type* = 0
      )
      {
        return return_helper<ix_t>(
          10 * ndtt_GC0<opts, dim>(psi, ndtt_GC[dim], i, j)
        );
      }
      
      template <opts_t opts, int dim, solvers::sptl_intrp_t sptl_intrp, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_spatial_helper(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<sptl_intrp == solvers::exact>::type* = 0
      )
      {
        return return_helper<ix_t>(
          ndxx_GC0<opts, dim>(psi, GC[dim], i, j)
        );
      }
      
      template <opts_t opts, int dim, solvers::sptl_intrp_t sptl_intrp, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_spatial_helper(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<sptl_intrp == solvers::aver2>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 * ndxx_GC0<opts, dim>(psi, GC[dim], i, j)
        );
      }
      
      template <opts_t opts, int dim, solvers::sptl_intrp_t sptl_intrp, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_spatial_helper(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<sptl_intrp == solvers::aver4>::type* = 0
      )
      {
        return 0;
      }
      
      template <opts_t opts, int dim, solvers::sptl_intrp_t sptl_intrp, class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd_spatial(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G, 
        const ix_t &i, 
        const ix_t &j
      )
      {
        return return_helper<ix_t>(
          - 1.0 / 24 *
          (
              4 * GC[dim](pi<dim>(i+h, j)) * ndxx_psi<opts, dim>(psi, i, j)
            + 2 * ndx_psi<opts, dim>(psi, i, j) * ndx_GC0<dim>(GC[dim], i, j)
            + div_3rd_spatial_helper<opts, dim, sptl_intrp>(psi, GC, i, j)
          )
        );
      }
      
      template <opts_t opts, int dim,
                solvers::sptl_intrp_t, solvers::tmprl_extrp_t,
                class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd(
        const arr_2d_t &psi_np1, 
        const arr_2d_t &psi_n, 
        const arrvec_t<arr_2d_t> &GC,
        const arrvec_t<arr_2d_t> &ndt_GC,
        const arrvec_t<arr_2d_t> &ndtt_GC,
        const arr_2d_t &G, 
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::div_3rd) && !opts::isset(opts, opts::div_3rd_dt)>::type* = 0
      )
      {
        return 0;
      }
      
      template <opts_t opts, int dim,
                solvers::sptl_intrp_t sptl_intrp, solvers::tmprl_extrp_t tmprl_extrp,
                class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd(
        const arr_2d_t &psi_np1, 
        const arr_2d_t &psi_n, 
        const arrvec_t<arr_2d_t> &GC,
        const arrvec_t<arr_2d_t> &ndt_GC,
        const arrvec_t<arr_2d_t> &ndtt_GC,
        const arr_2d_t &G, 
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::div_3rd)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          // upwind differencing correction
          div_3rd_upwind<opts, dim>(psi_np1, GC, G, i, j)
          // spatial terms
          + div_3rd_spatial<opts, dim, sptl_intrp>(psi_np1, GC, G, i, j)
          // mixed terms
          + 0.5 * abs(GC[dim](pi<dim>(i+h, j))) * ndx_fdiv<opts, dim>(psi_np1, GC, G, i, j)
          // temporal terms
          + 1.0 / 24 *
          (
              - 8 * GC[dim](pi<dim>(i+h, j)) *  nfdiv_fdiv<opts, dim>(psi_np1, GC, G, i, j)
              + div_3rd_temporal<opts, dim, tmprl_extrp>(psi_np1, ndtt_GC, i, j)
              + 2 * GC[dim](pi<dim>(i+h, j)) *  nfdiv<opts, dim>(psi_np1, ndt_GC, G, i, j)
              - 2 * ndt_GC[dim](pi<dim>(i+h, j)) * nfdiv<opts, dim>(psi_np1, GC, G, i, j)
          )
        );
      }
      
      template <opts_t opts, int dim,
                solvers::sptl_intrp_t sptl_intrp, solvers::tmprl_extrp_t tmprl_extrp,
                class arr_2d_t, class ix_t>
      forceinline_macro auto div_3rd(
        const arr_2d_t &psi_np1, 
        const arr_2d_t &psi_n, 
        const arrvec_t<arr_2d_t> &GC,
        const arrvec_t<arr_2d_t> &ndt_GC,
        const arrvec_t<arr_2d_t> &ndtt_GC,
        const arr_2d_t &G, 
        const ix_t &i, 
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::div_3rd_dt)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          // upwind differencing correction
          div_3rd_upwind<opts, dim>(psi_np1, GC, G, i, j)
          // spatial terms
          + div_3rd_spatial<opts, dim, sptl_intrp>(psi_np1, GC, G, i, j)
          // mixed terms
          - 0.5 * abs(GC[dim](pi<dim>(i+h, j))) * ndtx_psi<opts, dim>(psi_np1, psi_n, i, j)
          // temporal terms
          + 1.0 / 24 *
          (
              + 8 * GC[dim](pi<dim>(i+h, j)) *  nfdiv_dt<opts, dim>(psi_np1, psi_n, GC, G, i, j)
              + 1 * ndtt_GC0<opts, dim>(psi_np1, ndtt_GC[dim], i, j)
              + 2 * GC[dim](pi<dim>(i+h, j)) *  nfdiv<opts, dim>(psi_np1, ndt_GC, G, i, j)
              + 2 * ndt_GC[dim](pi<dim>(i+h, j)) * ndt_psi<opts, dim>(psi_np1, psi_n, i, j)
          )
        );
      }
      
      // antidiffusive velocity - standard version
      template <opts_t opts, int dim, solvers::sptl_intrp_t, solvers::tmprl_extrp_t, class arr_2d_t>
      inline void antidiff(
        arr_2d_t &res, 
        const arr_2d_t &psi_np1, 
        const arr_2d_t &psi_n, 
        const arrvec_t<arr_2d_t> &GC,
        const arrvec_t<arr_2d_t> &ndt_GC, // to have consistent interface with the div_3rd version
        const arrvec_t<arr_2d_t> &ndtt_GC, // ditto
        const arr_2d_t &G, 
        const rng_t &ir, 
        const rng_t &jr,
        typename std::enable_if<!opts::isset(opts, opts::div_2nd) && !opts::isset(opts, opts::div_3rd)>::type* = 0
      ) 
      {
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            res(pi<dim>(i, j)) = 
            // second order terms
            abs(GC[dim](pi<dim>(i+h, j))) / 2
            * (1 - abs(GC[dim](pi<dim>(i+h, j))) / G_bar_x<opts, dim>(G, i, j))
            * ndx_psi<opts, dim>(psi_np1, i, j) 
            - 
            GC[dim](pi<dim>(i+h, j)) 
            * GC1_bar_xy<dim>(GC[dim+1], i, j)
            / (2 * G_bar_x<opts, dim>(G, i, j))
            * ndy_psi<opts, dim>(psi_np1, i, j)
            // third order terms
            + TOT<opts, dim>(psi_np1, GC, G, i, j)
            //// fourth order terms
            + FOT<opts, dim>(psi_np1, GC, G, i, j)
            // divergent flow correction
            + DFL<opts, dim>(psi_np1, GC, G, i, j);
          }
        }
      }

      // antidiffusive velocity - divergence form
      template <opts_t opts, int dim, solvers::sptl_intrp_t sptl_intrp, solvers::tmprl_extrp_t tmprl_extrp, class arr_2d_t>
      inline void antidiff(
        arr_2d_t &res, 
        const arr_2d_t &psi_np1, 
        const arr_2d_t &psi_n, 
        const arrvec_t<arr_2d_t> &GC,
        const arrvec_t<arr_2d_t> &ndt_GC,
        const arrvec_t<arr_2d_t> &ndtt_GC,
        const arr_2d_t &G, 
        const rng_t &ir, 
        const rng_t &jr,
        typename std::enable_if<opts::isset(opts, opts::div_2nd)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::tot), "div_2nd & div_3rd options are incompatible with tot");
        static_assert(!opts::isset(opts, opts::dfl), "div_2nd & div_3rd options are incompatible with dfl");
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            res(pi<dim>(i + h, j)) = 
            div_2nd<opts, dim>(psi_np1, GC, G, i, j) +
            div_3rd<opts, dim, sptl_intrp, tmprl_extrp>(psi_np1, psi_n, GC, ndt_GC, ndtt_GC, G, i, j)
            // fourth order terms
            + FOT<opts, dim>(psi_np1, GC, G, i, j);
          }
        }
      } 
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx 
