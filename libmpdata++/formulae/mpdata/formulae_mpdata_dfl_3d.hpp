/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_3d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_psi_3d.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // divergent flow correction - no correction
      template<opts_t opts, int dim, class arr_3d_t>
      inline auto DFL(
        const arr_3d_t &psi,
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::dfl)>::type* = 0 
      ) -> decltype(0)
      { 
        return 0;  
      }

      // divergent flow correction - general case
      template<opts_t opts, int dim, class arr_3d_t>
      inline auto DFL(
        const arr_3d_t &psi,    //to have the same arguments as in iga option
        const arrvec_t<arr_3d_t> &GC,      
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::dfl) && !opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(,
        - 0.25 * GC[dim](pi<dim>(i+h, j, k)) 
        /
        G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        * 
        (
          (
            GC[dim](pi<dim>((i+1)+h, j, k)) - 
            GC[dim](pi<dim>(i-h    , j, k))
          )
          +
          (
            GC[dim+1](pi<dim>(i+1, j+h, k)) + 
            GC[dim+1](pi<dim>(i,   j+h, k)) -
            GC[dim+1](pi<dim>(i+1, j-h, k)) - 
            GC[dim+1](pi<dim>(i,   j-h, k))
          )
          +
          (
            GC[dim-1](pi<dim>(i+1, j, k+h)) + 
            GC[dim-1](pi<dim>(i,   j, k+h)) -
            GC[dim-1](pi<dim>(i+1, j, k-h)) - 
            GC[dim-1](pi<dim>(i,   j, k-h))
          )
        )
      )

      // divergent flow correction - infinite gauge version
      template<opts_t opts, int dim, class arr_3d_t>
      inline auto DFL(
        const arr_3d_t &psi,
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::dfl) && opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(,
        - 0.25 * GC[dim](pi<dim>(i+h, j, k)) 
        /
        G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        * 
        (
          (
            GC[dim](pi<dim>((i+1)+h, j, k)) - 
            GC[dim](pi<dim>(i-h    , j, k))
          )
          +
          (
            GC[dim+1](pi<dim>(i+1, j+h, k)) + 
            GC[dim+1](pi<dim>(i,   j+h, k)) -
            GC[dim+1](pi<dim>(i+1, j-h, k)) - 
            GC[dim+1](pi<dim>(i,   j-h, k))
          )
          +
          (
            GC[dim-1](pi<dim>(i+1, j, k+h)) + 
            GC[dim-1](pi<dim>(i,   j, k+h)) -
            GC[dim-1](pi<dim>(i+1, j, k-h)) - 
            GC[dim-1](pi<dim>(i,   j, k-h))
          )
        )
        * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
      )
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
