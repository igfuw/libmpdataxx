/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions that involve the flux divergence
// note that the function naming convention and the comments in this file correspond to dim = 0,
// if dim = 1 or 2 then you have to mentally perform appropiate substitutions
// for example for dim = 1 (x, y, z) -> (y, z, x) and (i+1/2, j, k) -> (i, j+1/2, k),
//             for dim = 2 (x, y, z) -> (z, x, y) and (i+1/2, j, k) -> (i, j, k+1/2)

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_psi_3d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_gc_3d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_3d.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // flux divergence i.e. 
      // 1 / G * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy + dz * d(GC * psi)/dz) at (i, j, k)
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto fdiv_centre(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            GC[dim  ](pi<dim>(i+h, j  , k  )) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i  , j  , k  )
          - GC[dim  ](pi<dim>(i-h, j  , k  )) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i-1, j  , k  )

          + GC[dim+1](pi<dim>(i  , j+h, k  )) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i  , j  , k  )
          - GC[dim+1](pi<dim>(i  , j-h, k  )) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i  , j-1, k  )

          + GC[dim-1](pi<dim>(i  , j  , k+h)) * psi_bar_z<opts BOOST_PP_COMMA() dim>(psi, i  , j  , k  )
          - GC[dim-1](pi<dim>(i  , j  , k-h)) * psi_bar_z<opts BOOST_PP_COMMA() dim>(psi, i  , j  , k-1)
          ) / formulae::G<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        );
      }
      
      // flux divergence i.e. 
      // 1 / G * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy + dz * d(GC * psi)/dz) at (i+1/2, j+1/2, k)
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto fdiv_corner_xy(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
              GC0_bar_xy<dim>(GC[dim  ], i+1, j, k) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i+1, j, k)
            - GC0_bar_xy<dim>(GC[dim  ], i  , j, k) * psi_bar_y<opts BOOST_PP_COMMA() dim>(psi, i  , j, k)
              
            + GC1_bar_xy<dim>(GC[dim+1], i, j+1, k) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j+1, k)
            - GC1_bar_xy<dim>(GC[dim+1], i, j  , k) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j  , k)
            
            + GC2_bar_xy<dim>(GC[dim-1], i, j, k  ) * psi_bar_xyz<opts BOOST_PP_COMMA() dim>(psi, i, j , k  )
            - GC2_bar_xy<dim>(GC[dim-1], i, j, k-1) * psi_bar_xyz<opts BOOST_PP_COMMA() dim>(psi, i, j , k-1)
          ) / G_bar_xy<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        );
      }
      
      // flux divergence i.e. 
      // 1 / G * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy + dz * d(GC * psi)/dz) at (i+1/2, j, k+1/2)
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto fdiv_corner_xz(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
              GC0_bar_xz<dim>(GC[dim  ], i+1, j, k) * psi_bar_z<opts BOOST_PP_COMMA() dim>(psi, i+1, j, k)
            - GC0_bar_xz<dim>(GC[dim  ], i  , j, k) * psi_bar_z<opts BOOST_PP_COMMA() dim>(psi, i  , j, k)
              
            + GC1_bar_xz<dim>(GC[dim+1], i, j  , k) * psi_bar_xyz<opts BOOST_PP_COMMA() dim>(psi, i, j  , k)
            - GC1_bar_xz<dim>(GC[dim+1], i, j-1, k) * psi_bar_xyz<opts BOOST_PP_COMMA() dim>(psi, i, j-1, k)
            
            + GC2_bar_xz<dim>(GC[dim-1], i, j, k+1) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j , k+1)
            - GC2_bar_xz<dim>(GC[dim-1], i, j, k  ) * psi_bar_x<opts BOOST_PP_COMMA() dim>(psi, i, j , k  )
          ) / G_bar_xz<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        );
      }
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy + dz * d(GC * psi)/dz) at (i+1/2, j, k)
      // positive sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          10 *
          frac<opts, ix_t>(
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) * psi(pi<dim>(i+1, j, k))
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) * psi(pi<dim>(i  , j, k))

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j  , k)
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j-1, k)
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi, i, j, k  )
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) * psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi, i, j, k-1)
          ,
            G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k) * (
              psi(pi<dim>(i+1, j  , k  ))
            + psi(pi<dim>(i  , j  , k  ))
            + psi(pi<dim>(i  , j+1, k  ))
            + psi(pi<dim>(i+1, j+1, k  ))
            + psi(pi<dim>(i  , j-1, k  ))
            + psi(pi<dim>(i+1, j-1, k  ))
            + psi(pi<dim>(i+1, j  , k+1))
            + psi(pi<dim>(i+1, j  , k-1))
            + psi(pi<dim>(i  , j  , k-1))
            + psi(pi<dim>(i  , j  , k+1))
            )
          )
        );
      }
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy + dz * d(GC * psi)/dz) at (i+1/2, j, k)
      // variable sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          10 *
          frac<opts, ix_t>(
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) * abs(psi(pi<dim>(i+1, j, k)))
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) * abs(psi(pi<dim>(i  , j, k)))

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j  , k)
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j-1, k)
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi, i, j, k  )
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) * psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi, i, j, k-1)
          ,
            G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k) * (
              abs(psi(pi<dim>(i+1, j  , k  )))
            + abs(psi(pi<dim>(i  , j  , k  )))
            + abs(psi(pi<dim>(i  , j+1, k  )))
            + abs(psi(pi<dim>(i+1, j+1, k  )))
            + abs(psi(pi<dim>(i  , j-1, k  )))
            + abs(psi(pi<dim>(i+1, j-1, k  )))
            + abs(psi(pi<dim>(i+1, j  , k+1)))
            + abs(psi(pi<dim>(i+1, j  , k-1)))
            + abs(psi(pi<dim>(i  , j  , k-1)))
            + abs(psi(pi<dim>(i  , j  , k+1)))
            )
          )
        );
      }
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * (dx * d(GC * psi)/dx + dy * d(GC * psi)/dy + dz * d(GC * psi)/dz) at (i+1/2, j, k)
      // infinite gauge version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) * psi(pi<dim>(i+1, j, k))
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) * psi(pi<dim>(i  , j, k))

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j  , k)
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) * psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi, i, j-1, k)
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi, i, j, k  )
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) * psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi, i, j, k-1)
          ) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        );
      }
      
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * (dx * d(GC * fdiv)/dx + dy * d(GC * fdiv)/dy + dz * d(GC * fdiv)/dz) at (i+1/2, j, k)
      // positive sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          12 *
          frac<opts, ix_t>(
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j, k)
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j, k)

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j  , k)
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j-1, k)
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * fdiv_corner_xz<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k  )
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) * fdiv_corner_xz<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k-1)
          ,
            G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k) * (
              psi(pi<dim>(i+2,   j, k)) +
              psi(pi<dim>(i+1,   j, k)) +
              psi(pi<dim>(i  ,   j, k)) +
              psi(pi<dim>(i-1,   j, k)) +
              psi(pi<dim>(i  , j+1, k)) +
              psi(pi<dim>(i  , j-1, k)) +
              psi(pi<dim>(i+1, j+1, k)) +
              psi(pi<dim>(i+1, j-1, k)) +
              psi(pi<dim>(i  , j, k+1)) +
              psi(pi<dim>(i  , j, k-1)) +
              psi(pi<dim>(i+1, j, k+1)) +
              psi(pi<dim>(i+1, j, k-1))
            )
          )
        );
      }
      
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * (dx * d(GC * fdiv)/dx + dy * d(GC * fdiv)/dy + dz * d(GC * fdiv)/dz) at (i+1/2, j, k)
      // variable sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          12 *
          frac<opts, ix_t>(
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j, k)
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j, k)

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j  , k)
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j-1, k)
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * fdiv_corner_xz<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k  )
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) * fdiv_corner_xz<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k-1)
          ,
            G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k) * (
              abs(psi(pi<dim>(i+2,   j, k))) +
              abs(psi(pi<dim>(i+1,   j, k))) +
              abs(psi(pi<dim>(i  ,   j, k))) +
              abs(psi(pi<dim>(i-1,   j, k))) +
              abs(psi(pi<dim>(i  , j+1, k))) +
              abs(psi(pi<dim>(i  , j-1, k))) +
              abs(psi(pi<dim>(i+1, j+1, k))) +
              abs(psi(pi<dim>(i+1, j-1, k))) +
              abs(psi(pi<dim>(i  , j, k+1))) +
              abs(psi(pi<dim>(i  , j, k-1))) +
              abs(psi(pi<dim>(i+1, j, k+1))) +
              abs(psi(pi<dim>(i+1, j, k-1)))
            )
          )
        );
      }
      
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * (dx * d(GC * fdiv)/dx + dy * d(GC * fdiv)/dy + dz * d(GC * fdiv)/dz) at (i+1/2, j, k)
      // infinite gauge version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j, k)
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) * fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j, k)

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j  , k)
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) * fdiv_corner_xy<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j-1, k)
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * fdiv_corner_xz<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k  )
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) * fdiv_corner_xz<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k-1)
          ) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        );
      }
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(fdiv)/dx at (i+1/2, j, k) - positive sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          12 * 
          frac<opts, ix_t>(
            fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j, k)
          - fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j, k)
          ,
            psi(pi<dim>(i+2,   j, k)) +
            psi(pi<dim>(i+1,   j, k)) +
            psi(pi<dim>(i  ,   j, k)) +
            psi(pi<dim>(i-1,   j, k)) +
            psi(pi<dim>(i  , j+1, k)) +
            psi(pi<dim>(i  , j-1, k)) +
            psi(pi<dim>(i  , j, k-1)) +
            psi(pi<dim>(i  , j, k+1)) +
            psi(pi<dim>(i+1, j+1, k)) +
            psi(pi<dim>(i+1, j-1, k)) +
            psi(pi<dim>(i+1, j, k-1)) +
            psi(pi<dim>(i+1, j, k+1))
          )
        );
      }
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(fdiv)/dx at (i+1/2, j, k) - variable sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          12 * 
          frac<opts, ix_t>(
            fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j, k)
          - fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j, k)
          ,
            abs(psi(pi<dim>(i+2,   j, k))) +
            abs(psi(pi<dim>(i+1,   j, k))) +
            abs(psi(pi<dim>(i  ,   j, k))) +
            abs(psi(pi<dim>(i-1,   j, k))) +
            abs(psi(pi<dim>(i  , j+1, k))) +
            abs(psi(pi<dim>(i  , j-1, k))) +
            abs(psi(pi<dim>(i  , j, k-1))) +
            abs(psi(pi<dim>(i  , j, k+1))) +
            abs(psi(pi<dim>(i+1, j+1, k))) +
            abs(psi(pi<dim>(i+1, j-1, k))) +
            abs(psi(pi<dim>(i+1, j, k-1))) +
            abs(psi(pi<dim>(i+1, j, k+1)))
          )
        );
      }
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(fdiv)/dx at (i+1/2, j, k) - infinite gauge version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_3d_t &psi,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i+1, j, k)
          - fdiv_centre<opts BOOST_PP_COMMA() dim>(psi, GC, G, i  , j, k)
        );
      }
      
      // nondimensionalised flux divergence of time derivative of psi i.e. 
      // 1 / (G * psi) * (dx * d(GC * dpsi/dt)/dx + dy * d(GC * dpsi/dt)/dy + dz * d(GC * dpsi/dt)/dz) at (i+1/2, j, k)
      // positive sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_dt(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          20 *
          frac<opts, ix_t>(
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) *
            (psi_np1(pi<dim>(i+1, j, k)) - psi_n(pi<dim>(i+1, j, k)))
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) *
            (psi_np1(pi<dim>(i  , j, k)) - psi_n(pi<dim>(i  , j, k)))

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) *
            (psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k) - psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k)) 
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) *
            (psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_np1, i, j-1, k) - psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_n, i, j-1, k)) 
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * 
            (psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k)- psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k))
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) *
            (psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k-1) - psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k-1))
          ,
            G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k) * (
              psi_np1(pi<dim>(i+1, j  , k  ))
            + psi_np1(pi<dim>(i  , j  , k  ))
            + psi_np1(pi<dim>(i  , j+1, k  ))
            + psi_np1(pi<dim>(i+1, j+1, k  ))
            + psi_np1(pi<dim>(i  , j-1, k  ))
            + psi_np1(pi<dim>(i+1, j-1, k  ))
            + psi_np1(pi<dim>(i+1, j  , k+1))
            + psi_np1(pi<dim>(i+1, j  , k-1))
            + psi_np1(pi<dim>(i  , j  , k-1))
            + psi_np1(pi<dim>(i  , j  , k+1))
            + psi_n(pi<dim>(i+1, j  , k  ))
            + psi_n(pi<dim>(i  , j  , k  ))
            + psi_n(pi<dim>(i  , j+1, k  ))
            + psi_n(pi<dim>(i+1, j+1, k  ))
            + psi_n(pi<dim>(i  , j-1, k  ))
            + psi_n(pi<dim>(i+1, j-1, k  ))
            + psi_n(pi<dim>(i+1, j  , k+1))
            + psi_n(pi<dim>(i+1, j  , k-1))
            + psi_n(pi<dim>(i  , j  , k-1))
            + psi_n(pi<dim>(i  , j  , k+1))
            )
          )
        );
      }

      // nondimensionalised flux divergence of time derivative of psi i.e. 
      // 1 / (G * psi) * (dx * d(GC * dpsi/dt)/dx + dy * d(GC * dpsi/dt)/dy + dz * d(GC * dpsi/dt)/dz) at (i+1/2, j, k)
      // variable sign scalar version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_dt(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          20 *
          frac<opts, ix_t>(
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) *
            (abs(psi_np1(pi<dim>(i+1, j, k))) - abs(psi_n(pi<dim>(i+1, j, k))))
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) *
            (abs(psi_np1(pi<dim>(i  , j, k))) - abs(psi_n(pi<dim>(i  , j, k))))

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) *
            (psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k) - psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k)) 
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) *
            (psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_np1, i, j-1, k) - psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_n, i, j-1, k)) 
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * 
            (psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k)- psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k))
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) *
            (psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k-1) - psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k-1))
          ,
            G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k) * (
              abs(psi_np1(pi<dim>(i+1, j  , k  )))
            + abs(psi_np1(pi<dim>(i  , j  , k  )))
            + abs(psi_np1(pi<dim>(i  , j+1, k  )))
            + abs(psi_np1(pi<dim>(i+1, j+1, k  )))
            + abs(psi_np1(pi<dim>(i  , j-1, k  )))
            + abs(psi_np1(pi<dim>(i+1, j-1, k  )))
            + abs(psi_np1(pi<dim>(i+1, j  , k+1)))
            + abs(psi_np1(pi<dim>(i+1, j  , k-1)))
            + abs(psi_np1(pi<dim>(i  , j  , k-1)))
            + abs(psi_np1(pi<dim>(i  , j  , k+1)))
            + abs(psi_n(pi<dim>(i+1, j  , k  )))
            + abs(psi_n(pi<dim>(i  , j  , k  )))
            + abs(psi_n(pi<dim>(i  , j+1, k  )))
            + abs(psi_n(pi<dim>(i+1, j+1, k  )))
            + abs(psi_n(pi<dim>(i  , j-1, k  )))
            + abs(psi_n(pi<dim>(i+1, j-1, k  )))
            + abs(psi_n(pi<dim>(i+1, j  , k+1)))
            + abs(psi_n(pi<dim>(i+1, j  , k-1)))
            + abs(psi_n(pi<dim>(i  , j  , k-1)))
            + abs(psi_n(pi<dim>(i  , j  , k+1)))
            )
          )
        );
      }
      
      // nondimensionalised flux divergence of time derivative of psi i.e. 
      // 1 / (G * psi) * (dx * d(GC * dpsi/dt)/dx + dy * d(GC * dpsi/dt)/dy + dz * d(GC * dpsi/dt)/dz) at (i+1/2, j, k)
      // infinite gauge version
      template <opts_t opts, int dim, class arr_3d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_dt(
        const arr_3d_t &psi_np1,
        const arr_3d_t &psi_n,
        const arrvec_t &GC, 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            GC0_bar_x<dim>(GC[dim  ], i+1, j, k) *
            (psi_np1(pi<dim>(i+1, j, k)) - psi_n(pi<dim>(i+1, j, k)))
          - GC0_bar_x<dim>(GC[dim  ], i  , j, k) *
            (psi_np1(pi<dim>(i  , j, k)) - psi_n(pi<dim>(i  , j, k)))

          + GC1_bar_x<dim>(GC[dim+1], i, j  , k) *
            (psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k) - psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k)) 
          - GC1_bar_x<dim>(GC[dim+1], i, j-1, k) *
            (psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_np1, i, j-1, k) - psi_bar_xy<opts BOOST_PP_COMMA() dim>(psi_n, i, j-1, k)) 
          
          + GC2_bar_x<dim>(GC[dim-1], i, j, k  ) * 
            (psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k)- psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k))
          - GC2_bar_x<dim>(GC[dim-1], i, j, k-1) *
            (psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_np1, i, j, k-1) - psi_bar_xz<opts BOOST_PP_COMMA() dim>(psi_n, i, j, k-1))
          ) / G_bar_x<opts BOOST_PP_COMMA() dim>(G, i, j, k)
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
