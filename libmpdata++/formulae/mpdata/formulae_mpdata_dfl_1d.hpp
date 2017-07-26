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

      //divergent flow correction see eq. (30) from @copybrief Smolarkiewicz_and_Margolin_1998)
      template<opts_t opts, class arr_1d_t, class ix_t>
      inline auto DFL(
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<!opts::isset(opts, opts::dfl)>::type* = 0 
      )
      { 
        return 0;  
      }

      template<opts_t opts, class arr_1d_t, class ix_t>
      inline auto DFL(
        const arr_1d_t &psi,    //to have the same arguments as in iga option
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::dfl) && !opts::isset(opts, opts::iga)>::type* = 0 
      )
      {
        return return_helper<ix_t>(
          - 0.5 * GC(i+h) 
          / 
          (formulae::G<opts>(G, i+1) + formulae::G<opts>(G, i)) 
          * 
          (GC((i+1)+h) - GC(i-h))
        );
      }

      template<opts_t opts, class arr_1d_t, class ix_t>
      inline auto DFL(
        const arr_1d_t &psi,
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const ix_t &i,
        typename std::enable_if<opts::isset(opts, opts::dfl) && opts::isset(opts, opts::iga)>::type* = 0 
      )
      {
        return return_helper<ix_t>(
          - 0.5 * GC(i+h) 
          / 
          (formulae::G<opts>(G, i+1) + formulae::G<opts>(G, i)) 
          * 
          (GC((i+1)+h) - GC(i-h))
          *
          0.5 *  (psi(i+1) + psi(i)) //to be compatible with iga formulation
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
