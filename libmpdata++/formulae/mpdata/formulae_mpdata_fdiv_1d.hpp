/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions that involve the flux divergence

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_psi_1d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_gc_1d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_1d.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // flux divergence i.e. 
      // 1 / G * dx * d(GC * psi)/dx at (i)
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto fdiv_centre(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i
      )
      {
        return return_helper<ix_t>(
          (
            GC[0](i+h) * psi_bar_x<opts>(psi, i  )
          - GC[0](i-h) * psi_bar_x<opts>(psi, i-1)
          ) / formulae::G<opts>(G, i)
        );
      }
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * dx * d(GC * psi)/dx at (i+1/2) - positive sign scalar version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto nfdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
            GC0_bar_x(GC[0], i+1) * psi(i+1)
          - GC0_bar_x(GC[0], i  ) * psi(i  )
          ,
            G_bar_x<opts>(G, i) * (
              psi(i+1)
            + psi(i  )
            )
          )
        );
      }
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * dx * d(GC * psi)/dx at (i+1/2) - variable sign scalar version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto nfdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          2 *
          frac<opts, ix_t>(
            GC0_bar_x(GC[0], i+1) * abs(psi(i+1))
          - GC0_bar_x(GC[0], i  ) * abs(psi(i  ))
          ,
            G_bar_x<opts>(G, i) * (
              abs(psi(i+1))
            + abs(psi(i  ))
            )
          )
        );
      }
      
      // nondimensionalised flux divergence i.e. 
      // 1 / (G * psi) * dx * d(GC * psi)/dx at (i+1/2) - infinite gauge version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto nfdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            GC0_bar_x(GC[0], i+1) * psi(i+1)
          - GC0_bar_x(GC[0], i  ) * psi(i  )
          ) / G_bar_x<opts>(G, i)
        );
      }
     
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * dx * d(GC * fdiv)/dx at (i+1/2) - positive sign scalar version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 *
          frac<opts, ix_t>(
            GC0_bar_x(GC[0], i+1) * fdiv_centre<opts>(psi, GC, G, i+1)
          - GC0_bar_x(GC[0], i  ) * fdiv_centre<opts>(psi, GC, G, i  )
          ,
            G_bar_x<opts>(G, i) * (
              psi(i+2) +
              psi(i+1) +
              psi(i  ) +
              psi(i-1)
            )
          )
        );
      }
      
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * dx * d(GC * fdiv)/dx at (i+1/2) - variable sign scalar version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 *
          frac<opts, ix_t>(
            GC0_bar_x(GC[0], i+1) * fdiv_centre<opts>(psi, GC, G, i+1)
          - GC0_bar_x(GC[0], i  ) * fdiv_centre<opts>(psi, GC, G, i  )
          ,
            G_bar_x<opts>(G, i) * (
              abs(psi(i+2)) +
              abs(psi(i+1)) +
              abs(psi(i  )) +
              abs(psi(i-1))
            )
          )
        );
      }
      
      // nondimensionalised divergence of flux divergence flux i.e.
      // 1 / (G * psi) * dx * d(GC * fdiv)/dx at (i+1/2) - infinite gauge version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto nfdiv_fdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          (
            GC0_bar_x(GC[0], i+1) * fdiv_centre<opts>(psi, GC, G, i+1)
          - GC0_bar_x(GC[0], i  ) * fdiv_centre<opts>(psi, GC, G, i  )
          ) / G_bar_x<opts>(G, i)
        );
      }
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(fdiv)/dx at (i+1/2) - positive sign scalar version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 * 
          frac<opts, ix_t>(
            fdiv_centre<opts>(psi, GC, G, i+1)
          - fdiv_centre<opts>(psi, GC, G, i  )
          ,
            psi(i+2) +
            psi(i+1) +
            psi(i  ) +
            psi(i-1)
          )
        );
      }
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(fdiv)/dx at (i+1/2) - variable sign scalar version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          4 * 
          frac<opts, ix_t>(
            fdiv_centre<opts>(psi, GC, G, i+1)
          - fdiv_centre<opts>(psi, GC, G, i  )
          ,
            abs(psi(i+2)) +
            abs(psi(i+1)) +
            abs(psi(i  )) +
            abs(psi(i-1))
          )
        );
      }
      
      // nondimensionalised x derivative of flux divergence i.e.
      // 1 / (G * psi) * dx * d(fdiv)/dx at (i+1/2) - infinite gauge version
      template <opts_t opts, class arr_1d_t, class ix_t, class arrvec_t>
      inline auto ndx_fdiv(
        const arr_1d_t &psi,
        const arrvec_t &GC, 
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return return_helper<ix_t>(
            fdiv_centre<opts>(psi, GC, G, i+1)
          - fdiv_centre<opts>(psi, GC, G, i  )
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
