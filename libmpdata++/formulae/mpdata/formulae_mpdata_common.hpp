/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/idxperm.hpp>
#include <libmpdata++/arakawa_c.hpp>

namespace libmpdataxx
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      using namespace arakawa_c;
      using idxperm::pi;

      const int halo = 1, n_tlev = 2;
// TODO
/*
  struct varsgn
  {
    template<class T> auto aon(const T &x) -> decltype(abs(x))
    {
      return abs(x);
    }
  };

  struct posdef
  {
    template<class T> T aon(const T &x)
    {
      return x;
    }
  };
*/

      template<class nom_t, class den_t>
      inline auto frac(
        const nom_t &nom, const den_t &den
      ) return_macro(,
        where(den > 0, nom / den, 0)
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
