/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

// various numerical expressions relating to the G factor
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
      // interpolation of G to (i+1/2, j, k) - general case
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto G_bar_x(
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
      )
      {
        return return_helper<ix_t>(
          (
            formulae::G<opts, dim>(G, i+1, j, k)
          + formulae::G<opts, dim>(G, i  , j, k)
          ) / 2
        );
      }

      // interpolation of G to (i+1/2, j, k) - constant G version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto G_bar_x( 
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
      )
      {
          return 1;
      }
      
      // interpolation of G to (i+1/2, j+1/2, k) - general case
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto G_bar_xy(
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
      )
      {
        return return_helper<ix_t>(
          (
            formulae::G<opts, dim>(G, i  , j  , k) +
            formulae::G<opts, dim>(G, i  , j+1, k) +
            formulae::G<opts, dim>(G, i+1, j  , k) +
            formulae::G<opts, dim>(G, i+1, j+1, k)
          ) / 4
        );
      }
      
      // interpolation of G to (i+1/2, j+1/2, k) - constant G version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto G_bar_xy(
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
      )
      {
        return 1;
      }
      
      // interpolation of G to (i+1/2, j, k+1/2) - general case
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto G_bar_xz(
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
      )
      {
        return return_helper<ix_t>(
          (
            formulae::G<opts, dim>(G, i  , j, k  ) +
            formulae::G<opts, dim>(G, i  , j, k+1) +
            formulae::G<opts, dim>(G, i+1, j, k  ) +
            formulae::G<opts, dim>(G, i+1, j, k+1)
          ) / 4
        );
      }
      
      // interpolation of G to (i+1/2, j, k+1/2) - constant G version
      template<opts_t opts, int dim, class arr_3d_t, class ix_t>
      inline auto G_bar_xz(
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k,
        typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
      ) {
        return 1;
      }
    } // namespace mpdata
  } // namespace formulae
} // namespace libmpdataxx 
