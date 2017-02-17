/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
//#include <libmpdata++/formulae/mpdata/formulae_mpdata_dfl_3d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_hot_3d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fdiv_3d.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      // first come helpers for divergence form of antidiffusive velocity
      template <opts_t opts, int dim, class arr_3d_t>
      inline auto div_2nd(
        const arr_3d_t &psi, 
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G, 
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
        // second order terms
        abs(GC[dim](pi<dim>(i+h, j, k))) / 2
        * ndx_psi<opts BOOST_PP_COMMA() dim>(psi, i, j, k) 
        - 
        GC[dim](pi<dim>(i+h, j, k)) / 2
        * nfdiv<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)
      )
      
      template <opts_t opts, int dim, class arr_3d_t>
      inline auto div_3rd_helper(
        const arr_3d_t &psi, 
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G, 
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
        abs(div_2nd<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)) / 2
        * ndx_psi<opts BOOST_PP_COMMA() dim>(psi, i, j, k) 
      )

      template <opts_t opts, int dim, class arr_3d_t>
      inline typename arr_3d_t::T_numtype div_3rd_helper(
        const arr_3d_t &psi, 
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G, 
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return 0;
      }
      
      template <opts_t opts, int dim, class arr_3d_t>
      inline typename arr_3d_t::T_numtype div_3rd(
        const arr_3d_t &psi, 
        const arrvec_t<arr_3d_t> &GC,
        const arrvec_t<arr_3d_t> &ndt_GC,
        const arrvec_t<arr_3d_t> &ndtt_GC,
        const arr_3d_t &G, 
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::div_3rd)>::type* = 0
      )
      {
        return 0;
      }
      
      template <opts_t opts, int dim, class arr_3d_t>
      inline auto div_3rd(
        const arr_3d_t &psi, 
        const arrvec_t<arr_3d_t> &GC,
        const arrvec_t<arr_3d_t> &ndt_GC,
        const arrvec_t<arr_3d_t> &ndtt_GC,
        const arr_3d_t &G, 
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::div_3rd)>::type* = 0
      ) return_macro(,
        // upwind differencing correction
        div_3rd_helper<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)
        // spatial terms
        - 1.0 / 24 *
        (
            4 * GC[dim](pi<dim>(i+h, j, k)) * ndxx_psi<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
          + 2 * ndx_psi<opts BOOST_PP_COMMA() dim>(psi, i, j, k) * ndx_GC0<dim>(GC[dim], i, j, k)
          + 1 * ndxx_GC0<opts BOOST_PP_COMMA() dim>(psi, GC[dim], i, j, k)
        )
        // mixed terms
        + 0.5 * abs(GC[dim](pi<dim>(i+h, j, k))) * ndx_fdiv<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)
        // temporal terms
        + 1.0 / 24 *
        (
            - 8 * GC[dim](pi<dim>(i+h, j, k)) *  nfdiv_fdiv<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)
            + 1 * ndtt_GC0<opts BOOST_PP_COMMA() dim>(psi, ndtt_GC[dim], i, j, k)
            + 2 * GC[dim](pi<dim>(i+h, j, k)) *  nfdiv<opts BOOST_PP_COMMA() dim>(psi, ndt_GC, G, i, j, k)
            - 2 * ndt_GC[dim](pi<dim>(i+h, j, k)) * nfdiv<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)
        )
      )

      // antidiffusive velocity - standard version
      template <opts_t opts, int dim, class arr_3d_t>
      inline auto antidiff(
        const arr_3d_t &psi,
        const arrvec_t<arr_3d_t> &GC,
        const arrvec_t<arr_3d_t> &ndt_GC, // to have consistent interface with the div_3rd version
        const arrvec_t<arr_3d_t> &ndtt_GC, // ditto
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::div_2nd) && !opts::isset(opts, opts::div_3rd)>::type* = 0
      ) return_macro(,
          // second order terms
          abs(GC[dim](pi<dim>(i+h, j, k))) / 2
        * (1 - abs(GC[dim](pi<dim>(i+h, j, k))) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k))
        * ndx_psi<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
        - GC[dim](pi<dim>(i+h, j, k)) / 2
        * (
            GC1_bar_xy<dim>(GC[dim+1], i, j, k)
          * ndy_psi<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
          + GC2_bar_xz<dim>(GC[dim-1], i, j, k)
          * ndz_psi<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
          )
          / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k)
          // third order terms
        + TOT<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)
        // divergent flow correction
        //+ DFL<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k) 
      )

      // antidiffusive velocity - divergence form
      template <opts_t opts, int dim, class arr_3d_t>
      inline auto antidiff(
        const arr_3d_t &psi, 
        const arrvec_t<arr_3d_t> &GC,
        const arrvec_t<arr_3d_t> &ndt_GC,
        const arrvec_t<arr_3d_t> &ndtt_GC,
        const arr_3d_t &G, 
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::div_2nd)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::tot), "div_2nd & div_3rd options are incompatible with tot");
        static_assert(!opts::isset(opts, opts::dfl), "div_2nd & div_3rd options are incompatible with dfl");
        ,
        div_2nd<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k) +
        div_3rd<opts BOOST_PP_COMMA() dim>(psi, GC, ndt_GC, ndtt_GC, G, i, j, k)
      ) 
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
