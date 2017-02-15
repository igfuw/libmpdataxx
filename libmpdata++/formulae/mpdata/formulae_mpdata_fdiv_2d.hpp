/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions that involve the flux divergence
// note that the function naming convention and comments in this file correspond to dim = 0,
// if dim = 1 then you have to mentally perform appropiate substitutions
// for example (x, y) -> (y, x) and (i+1/2, j) -> (i, j+1/2)

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_psi_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_gc_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_2d.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // flux divergence i.e. 
      // 1 / G * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy) at (i, j)
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto fdiv_centre(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
        (
          GC[dim](pi<dim>(i+h,   j)) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i  , j  )
        - GC[dim](pi<dim>(i-h,   j)) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i-1, j  )
        + GC[dim+1](pi<dim>(i, j+h)) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i  , j  )
        - GC[dim+1](pi<dim>(i, j-h)) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i  , j-1)
        ) / formulae::G<opts BOOST_PP_COMMA() dim>(G, i, j)
      )
      
      // flux divergence i.e. 
      // 1 / G * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy) at (i+1/2, j+1/2)
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto fdiv_corner_xy(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
        (
            GC0_bar_xy<dim>(GC[dim], i+1, j  ) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i+1, j)
          - GC0_bar_xy<dim>(GC[dim], i  , j  ) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i  , j)
            
          + GC1_bar_xy<dim>(GC[dim+1], i, j+1) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j+1)
          - GC1_bar_xy<dim>(GC[dim+1], i, j  ) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j  )
        ) / G_bar_xy<opts BOOST_PP_COMMA() dim>(G, i, j)
      )
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy) at (i+1/2, j) - positive sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto nfdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        6 *
        frac<opts>(
          GC0_bar_x<dim>(GC[dim]  , i+1, j  ) * psi(pi<dim>(i+1, j))
        - GC0_bar_x<dim>(GC[dim]  , i  , j  ) * psi(pi<dim>(i  , j))
        + GC1_bar_x<dim>(GC[dim+1], i  , j  ) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j  )
        - GC1_bar_x<dim>(GC[dim+1], i  , j-1) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j-1)
        ,
          G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j) * (
            psi(pi<dim>(i+1, j  ))
          + psi(pi<dim>(i  , j  ))
          + psi(pi<dim>(i  , j+1))
          + psi(pi<dim>(i+1, j+1))
          + psi(pi<dim>(i  , j-1))
          + psi(pi<dim>(i+1, j-1))
          )
        )
      )
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy) at (i+1/2, j) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto nfdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        6 *
        frac<opts>(
          GC0_bar_x<dim>(GC[dim]  , i+1, j  ) * abs(psi(pi<dim>(i+1, j)))
        - GC0_bar_x<dim>(GC[dim]  , i  , j  ) * abs(psi(pi<dim>(i  , j)))
        + GC1_bar_x<dim>(GC[dim+1], i  , j  ) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j  )
        - GC1_bar_x<dim>(GC[dim+1], i  , j-1) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j-1)
        ,
          G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j) * (
            abs(psi(pi<dim>(i+1, j  )))
          + abs(psi(pi<dim>(i  , j  )))
          + abs(psi(pi<dim>(i  , j+1)))
          + abs(psi(pi<dim>(i+1, j+1)))
          + abs(psi(pi<dim>(i  , j-1)))
          + abs(psi(pi<dim>(i+1, j-1)))
          )
        )
      )
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy) at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto nfdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
        (
          GC0_bar_x<dim>(GC[dim]  , i+1, j  ) * psi(pi<dim>(i+1, j))
        - GC0_bar_x<dim>(GC[dim]  , i  , j  ) * psi(pi<dim>(i  , j))
        + GC1_bar_x<dim>(GC[dim+1], i  , j  ) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j  )
        - GC1_bar_x<dim>(GC[dim+1], i  , j-1) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j-1)
        ) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j)
      )
     
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * (dx * d(GC * fdiv)/dx + dy * d(GC * fdiv)/dy) at (i+1/2, j) - positive sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        8 *
        frac<opts>(
          GC0_bar_x<dim>(GC[dim]  , i+1, j  ) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j)
        - GC0_bar_x<dim>(GC[dim]  , i  , j  ) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j)
        + GC1_bar_x<dim>(GC[dim+1], i  , j  ) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j  )
        - GC1_bar_x<dim>(GC[dim+1], i  , j-1) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j-1)
        ,
          G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j) * (
            psi(pi<dim>(i+2,   j)) +
            psi(pi<dim>(i+1,   j)) +
            psi(pi<dim>(i  ,   j)) +
            psi(pi<dim>(i-1,   j)) +
            psi(pi<dim>(i  , j+1)) +
            psi(pi<dim>(i  , j-1)) +
            psi(pi<dim>(i+1, j+1)) +
            psi(pi<dim>(i+1, j-1))
          )
        )
      )
      
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * (dx * d(GC * fdiv)/dx + dy * d(GC * fdiv)/dy) at (i+1/2, j) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        8 *
        frac<opts>(
          GC0_bar_x<dim>(GC[dim]  , i+1, j  ) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j)
        - GC0_bar_x<dim>(GC[dim]  , i  , j  ) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j)
        + GC1_bar_x<dim>(GC[dim+1], i  , j  ) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j  )
        - GC1_bar_x<dim>(GC[dim+1], i  , j-1) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j-1)
        ,
          G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j) * (
            abs(psi(pi<dim>(i+2,   j))) +
            abs(psi(pi<dim>(i+1,   j))) +
            abs(psi(pi<dim>(i  ,   j))) +
            abs(psi(pi<dim>(i-1,   j))) +
            abs(psi(pi<dim>(i  , j+1))) +
            abs(psi(pi<dim>(i  , j-1))) +
            abs(psi(pi<dim>(i+1, j+1))) +
            abs(psi(pi<dim>(i+1, j-1)))
          )
        )
      )
      
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * (dx * d(GC * fdiv)/dx + dy * d(GC * fdiv)/dy) at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
        (
          GC0_bar_x<dim>(GC[dim]  , i+1, j  ) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j)
        - GC0_bar_x<dim>(GC[dim]  , i  , j  ) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j)
        + GC1_bar_x<dim>(GC[dim+1], i  , j  ) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j  )
        - GC1_bar_x<dim>(GC[dim+1], i  , j-1) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j-1)
        ) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j)
      )
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(GC * fdiv)/dx at (i+1/2, j) - positive sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        8 * 
        frac<opts>(
          fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j)
        - fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j)
        ,
          psi(pi<dim>(i+2,   j)) +
          psi(pi<dim>(i+1,   j)) +
          psi(pi<dim>(i  ,   j)) +
          psi(pi<dim>(i-1,   j)) +
          psi(pi<dim>(i  , j+1)) +
          psi(pi<dim>(i  , j-1)) +
          psi(pi<dim>(i+1, j+1)) +
          psi(pi<dim>(i+1, j-1))
        )
      )
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(GC * fdiv)/dx at (i+1/2, j) - variable sign scalar version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        8 * 
        frac<opts>(
          fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j)
        - fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j)
        ,
          abs(psi(pi<dim>(i+2, j  ))) +
          abs(psi(pi<dim>(i+1, j  ))) +
          abs(psi(pi<dim>(i  , j  ))) +
          abs(psi(pi<dim>(i-1, j  ))) +
          abs(psi(pi<dim>(i  , j+1))) +
          abs(psi(pi<dim>(i  , j-1))) +
          abs(psi(pi<dim>(i+1, j+1))) +
          abs(psi(pi<dim>(i+1, j-1)))
        )
      )
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(GC * fdiv)/dx at (i+1/2, j) - infinite gauge version
      template <opts_t opts, int dim, class arr_2d_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_2d_t &psi,
        const arrvec_t &GC, 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
          fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j)
        - fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j)
      )
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
