/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <advoocat/blitz.hpp>
#include <advoocat/bcond/bcond.hpp>

namespace advoocat
{
  namespace bcond
  {
    template <typename real_t>
    class cyclic_left_common : public bcond_t<real_t>
    {
      protected:

      // member fields
      rng_t left_halo, rght_edge;

      public:

      // ctor
      cyclic_left_common(const rng_t &i, int halo) :
	left_halo(i.first() - halo    , i.first() - 1       ),
	rght_edge(i.last()  - halo + 1, i.last()            )
      {} 
    };

    template <typename real_t>
    class cyclic_rght_common : public bcond_t<real_t>
    {
      protected:

      // member fields
      rng_t left_edge, rght_halo;

      public:

      // ctor
      cyclic_rght_common(const rng_t &i, int halo) :
	rght_halo(i.last()  + 1       , i.last()  + halo    ),
	left_edge(i.first()           , i.first() + halo - 1)
      {} 
    };

  }; // namespace bcond
}; // namespace advoocat
