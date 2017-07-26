/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_psi_1d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_g_1d.hpp>

namespace libmpdataxx
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      template<opts_t opts, class arr_1d_t, class ix_t>
      inline auto ndxx_psi_coeff(
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const ix_t &i
      )
      {
        return return_helper<ix_t>(
          (
            3 * GC(i+h) * abs(GC(i+h)) / G_bar_x<opts>(G, i)
            - 2 * pow(GC(i+h), 3) / pow(G_bar_x<opts>(G, i), 2)  
            - GC(i+h)
          ) / 6
        );
      }

      // third order terms
      template<opts_t opts, class arr_1d_t, class ix_t>
      inline auto TOT(
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::tot)>::type* = 0 
      )
      {
        return return_helper<ix_t>(
            ndxx_psi<opts>(psi, i) * ndxx_psi_coeff<opts>(GC, G, i)
        );
      }
      
      template<opts_t opts, class arr_1d_t, class ix_t>
      inline auto TOT(
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::tot)>::type* = 0 
      )
      { 
        return 0;
      }
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
