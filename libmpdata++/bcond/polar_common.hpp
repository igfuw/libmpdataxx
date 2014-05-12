/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/bcond/bcond.hpp>
#include <libmpdata++/formulae/arakawa_c.hpp>

// TODO: move to detail

namespace libmpdataxx
{
  namespace bcond
  {
    using namespace arakawa_c;

// TODO: add detail namespace, move to detail directory?
    template <typename real_t>
    class polar_common : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using parent_t::parent_t; // inheriting ctor

      protected:

      // member fields
      const int pole;

      int polar_neighbours(const int j)
      {
        return (j + pole) % (2 * pole);
      }

      public:

      // ctor
      polar_common(const rng_t &i, const int halo, const int pole) :
        parent_t(i, halo),
        pole(pole)
      {} 
    };
  }; // namespace bcond
}; // namespace libmpdataxx
