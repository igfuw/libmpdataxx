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
      template<opts_t opts, int dim, class arr_2d_t>
      inline auto ndxx_psi_coeff(
        const arr_2d_t &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
        (
          3 * GC(pi<dim>(i+h, j)) * abs(GC(pi<dim>(i+h, j))) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j)
          - 2 * pow(GC(pi<dim>(i+h, j)), 3) / pow(G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j), 2)
          - GC(pi<dim>(i+h, j))
        ) / 6
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto ndxy_psi_coeff(
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
        (abs(GC[dim](pi<dim>(i+h, j))) - 2 * pow(GC[dim](pi<dim>(i+h, j)), 2) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j)) 
         * GC1_bar_xy<dim>(GC[dim+1], i, j) / (2 * G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j))
      )
      
      // third order terms
      template<opts_t opts, int dim, class arr_2d_t>
      inline auto TOT(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::tot)>::type* = 0
      ) return_macro(,
         ndxx_psi<opts BOOST_PP_COMMA() dim>(psi, i, j) * ndxx_psi_coeff<opts BOOST_PP_COMMA() dim>(GC[dim], G, i, j) 
         + 
         ndxy_psi<opts BOOST_PP_COMMA() dim>(psi, i, j) * ndxy_psi_coeff<opts BOOST_PP_COMMA() dim>(GC, G, i, j)
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline typename arr_2d_t::T_numtype TOT(
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::tot)>::type* = 0
      )
      {
        return 0;
      }
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx 
