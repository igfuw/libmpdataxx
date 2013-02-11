/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "../blitz.hpp"

namespace advoocat
{
  namespace bcond
  {
    template <int n_dims, typename real_t = float>
    class cyclic_common
    {
      public:

      typedef blitz::Array<real_t, n_dims> arr_t;

      protected:

      // member fields
      rng_t left_halo, rght_halo;
      rng_t left_edge, rght_edge;

      public:

      // ctor
      cyclic_common(const rng_t &i, int halo) :
	left_halo(i.first() - halo    , i.first() - 1       ),
	rght_edge(i.last()  - halo + 1, i.last()            ),
	rght_halo(i.last()  + 1       , i.last()  + halo    ),
	left_edge(i.first()           , i.first() + halo - 1)
      {} 
    };
  }; // namespace bcond
}; // namespace advoocat
