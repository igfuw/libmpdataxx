/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/formulae/idxperm.hpp>
#include <libmpdata++/formulae/arakawa_c.hpp>
#include <libmpdata++/formulae/common.hpp>

//#include <boost/preprocessor/control/if.hpp>

namespace libmpdataxx
{ 
  // in namespace solvers to be consistent with other options
  namespace solvers
  {
    enum sptl_intrp_t
    { 
      exact,
      aver2,
      aver4,
    };
    
    enum tmprl_extrp_t
    { 
      noextrp,
      linear2,
    };
  }

  namespace formulae 
  { 
    namespace mpdata 
    {
      using namespace arakawa_c;
      using idxperm::pi;
      using opts::opts_t;
      using std::abs;

      const int n_tlev = 2;

      constexpr const int halo(const opts_t &opts) 
      {
        return (
          opts::isset(opts, opts::tot)        || // see psi 2-nd derivatives in eq. (36) in PKS & LGM 1998
          opts::isset(opts, opts::dfl)        || // see +3/2 in eq. (30) in PKS & LGM 1998
          opts::isset(opts, opts::div_2nd)    || 
          opts::isset(opts, opts::div_3rd)    ||
          opts::isset(opts, opts::div_3rd_dt)
        ) ? 2 : 1; 
      }

      // frac: implemented using blitz::where()
      template<opts_t opts, class ix_t, class nom_t, class den_t>
      forceinline_macro auto frac(
        const nom_t &nom, 
        const den_t &den,
        typename std::enable_if<opts::isset(opts, opts::pfc)>::type* = 0 // enabled if pfc == true
      )
      {
        return return_helper<ix_t>(
          where(den != 0, nom / den, 0) // note: apparently >0 would suffice
        );
      }

      // frac: implemented as suggested in MPDATA papers
      //       if den == 0, then adding a smallest representable positive number
      template<opts_t opts, class ix_t, class nom_t, class den_t>
      forceinline_macro auto frac(
        const nom_t &nom, 
        const den_t &den,
        typename std::enable_if<!opts::isset(opts, opts::pfc)>::type* = 0 // enabled if pfc == false
      )
      {
        return return_helper<ix_t>(
          nom / (den + blitz::tiny(typename real_t_helper<ix_t, nom_t>::type(0.))) // note: for negative signal eps -> -eps
        );
      }

      // a bigger-epsilon version for FCT (used regardless of opts::eps setting)
      template<class ix_t, class nom_t, class den_t>
      forceinline_macro auto fct_frac(
        const nom_t &nom, 
        const den_t &den
      )
      {
        return return_helper<ix_t>(
          nom / (den + blitz::epsilon(typename real_t_helper<ix_t, nom_t>::type(0.)))
        );
      }
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
