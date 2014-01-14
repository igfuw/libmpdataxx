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
      template<opts_t opts, int dim, class arr_2d_t>
      inline typename arr_2d_t::T_numtype G_at_half(  
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == false
      ) {
        return 1;
      }

      //G evaluated at (1+1/2, j)
      template<opts_t opts, int dim, class arr_2d_t>
      inline auto G_at_half( 
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::nug)>::type* = 0 // enabled if nug == true
      ) return_macro(,
        (
          formulae::G<opts>(G, pi<dim>(i+1, j)) 
          +
          formulae::G<opts>(G, pi<dim>(i  , j))
        ) / 2
      )

      //eq. (29b) from @copybrief Smolarkiewicz_and_Margolin_1998
      template<int dim, class arr_2d_t>
      inline auto GC_bar( // caution proper call looks like C_bar<dim>(C[dim-1], i, j) - note dim vs dim-1
        const arr_2d_t &GC, 
        const rng_t &i, 
        const rng_t &j
      ) return_macro(,
        (
          GC(pi<dim>(i+1, j+h)) + 
          GC(pi<dim>(i,   j+h)) +
          GC(pi<dim>(i+1, j-h)) + 
          GC(pi<dim>(i,   j-h)) 
        ) / 4
      )

    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx 
