/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_dfl_1d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_hot_1d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_fdiv_1d.hpp>

namespace libmpdataxx
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // first come helpers for divergence form of antidiffusive velocity
      template <opts_t opts, class arr_1d_t>
      inline auto div_2nd(
        const arr_1d_t &psi, 
        const arrvec_t<arr_1d_t> &GC,
        const arr_1d_t &G, 
        const rng_t &i
      ) return_macro(,
        abs(GC[0](i+h)) / 2
        * ndx_psi<opts>(psi, i) 
        - 
        GC[0](i+h) / 2
        * nfdiv<opts>(psi, GC, G, i)
      )
      
      template <opts_t opts, class arr_1d_t>
      inline auto div_3rd_helper(
        const arr_1d_t &psi, 
        const arrvec_t<arr_1d_t> &GC,
        const arr_1d_t &G, 
        const rng_t &i, 
        typename std::enable_if<!opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
        abs(div_2nd<opts>(psi, GC, G, i)) / 2
        * ndx_psi<opts>(psi, i) 
      )

      template <opts_t opts, class arr_1d_t>
      inline typename arr_1d_t::T_numtype div_3rd_helper(
        const arr_1d_t &psi, 
        const arrvec_t<arr_1d_t> &GC,
        const arr_1d_t &G, 
        const rng_t &i, 
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      )
      {
        return 0;
      }
      
      template <opts_t opts, class arr_1d_t>
      inline typename arr_1d_t::T_numtype div_3rd(
        const arr_1d_t &psi, 
        const arrvec_t<arr_1d_t> &GC,
        const arrvec_t<arr_1d_t> &ndt_GC,
        const arrvec_t<arr_1d_t> &ndtt_GC,
        const arr_1d_t &G, 
        const rng_t &i, 
        typename std::enable_if<!opts::isset(opts, opts::div_3rd)>::type* = 0
      )
      {
        return 0;
      }
      
      template <opts_t opts, class arr_1d_t>
      inline auto div_3rd(
        const arr_1d_t &psi, 
        const arrvec_t<arr_1d_t> &GC,
        const arrvec_t<arr_1d_t> &ndt_GC,
        const arrvec_t<arr_1d_t> &ndtt_GC,
        const arr_1d_t &G, 
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::div_3rd)>::type* = 0
      ) return_macro(,
        // upwind differencing correction
        div_3rd_helper<opts>(psi, GC, G, i)
        // spatial terms
        - 1.0 / 24 *
        (
            4 * GC[0](i+h) * ndxx_psi<opts>(psi, i)
          + 2 * ndx_psi<opts>(psi, i) * ndx_GC0(GC[0], i)
          + 1 * ndxx_GC0<opts>(psi, GC[0], i)
        )
        // mixed terms
        + 0.5 * abs(GC[0](i+h)) * ndx_fdiv<opts>(psi, GC, G, i)
        // temporal terms
        + 1.0 / 24 *
        (
            - 8 * GC[0](i+h) *  nfdiv_fdiv<opts>(psi, GC, G, i)
            + 1 * ndtt_GC0<opts>(psi, ndtt_GC[0], i)
            + 2 * GC[0](i+h) *  nfdiv<opts>(psi, ndt_GC, G, i)
            - 2 * ndt_GC[0](i+h) * nfdiv<opts>(psi, GC, G, i)
        )
      )
      
      // antidiffusive velocity - standard version
      template<opts_t opts, class arr_1d_t>
      inline void antidiff( // antidiffusive velocity
        arr_1d_t &res, 
        const arr_1d_t &psi, 
        const arrvec_t<arr_1d_t> &GC,
        const arrvec_t<arr_1d_t> &ndt_GC, // to have consistent interface with the div_3rd version
        const arrvec_t<arr_1d_t> &ndtt_GC, // ditto
        const arr_1d_t &G,
        const rng_t &ir,
        typename std::enable_if<!opts::isset(opts, opts::div_2nd) && !opts::isset(opts, opts::div_3rd)>::type* = 0
      )
      {
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          res(i) = 
          // second-order terms
          abs(GC[0](i+h)) / 2
          * (1 - abs(GC[0](i+h)) / G_bar_x<opts>(G, i))
          * ndx_psi<opts>(psi, i) 
          // third-order terms
          + TOT<opts>(psi, GC[0], G, i) //higher order term
          // divergent flow terms
          + DFL<opts>(psi, GC[0], G, i); //divergent flow correction
        }
      }

      // antidiffusive velocity - divergence form
      template <opts_t opts, class arr_1d_t>
      inline void antidiff(
        arr_1d_t &res,
        const arr_1d_t &psi, 
        const arrvec_t<arr_1d_t> &GC,
        const arrvec_t<arr_1d_t> &ndt_GC,
        const arrvec_t<arr_1d_t> &ndtt_GC,
        const arr_1d_t &G, 
        const rng_t &ir, 
        typename std::enable_if<opts::isset(opts, opts::div_2nd)>::type* = 0
      )
      {
        static_assert(!opts::isset(opts, opts::tot), "div_2nd & div_3rd options are incompatible with tot");
        static_assert(!opts::isset(opts, opts::dfl), "div_2nd & div_3rd options are incompatible with dfl");

        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          res(i) = 
          div_2nd<opts>(psi, GC, G, i) +
          div_3rd<opts>(psi, GC, ndt_GC, ndtt_GC, G, i);
        }
      } 

    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
