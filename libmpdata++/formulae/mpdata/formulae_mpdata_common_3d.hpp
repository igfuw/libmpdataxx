/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      template<opts_t opts, int dim, class arr_3d_t>
      inline typename arr_3d_t::T_numtype G_at_half( 
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
      ) {
          return 1;
      }

      //G evaluated at (i+1/2, j, k)
      template<opts_t opts, int dim, class arr_3d_t>
      inline auto G_at_half(
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
      ) return_macro(,
          (
            formulae::G<opts, dim>(G, i+1, j, k)
          + formulae::G<opts, dim>(G, i  , j, k)
          ) / 2
      )

      template<int d, class arr_3d_t>
      inline auto GC_bar1( // caution proper call looks like GC_bar1<dim>(GC[dim+1], i, j, k) - note dim vs dim+1
        const arr_3d_t &GC,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
          (
            GC(pi<d>(i+1, j+h, k))
          + GC(pi<d>(i,   j+h, k))
          + GC(pi<d>(i+1, j-h, k))
          + GC(pi<d>(i,   j-h, k))
          ) / 4
      )
     
      template<int d, class arr_3d_t>
      inline auto GC_bar2( // caution proper call looks like GC_bar2<dim>(GC[dim-1], i, j, k) - note dim vs dim-1
        const arr_3d_t &GC,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
          (
            GC(pi<d>(i+1, j, k+h))
          + GC(pi<d>(i,   j, k+h))
          + GC(pi<d>(i+1, j, k+h))
          + GC(pi<d>(i,   j, k+h))
          ) / 4
      )

    }; // namespace mpdata
  }; // namespace formulae
}; // namespace libmpdataxx
