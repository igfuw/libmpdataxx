/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/bcond/bcond.hpp>
#include <libmpdata++/arakawa_c.hpp>

// TODO: move to detail

namespace libmpdataxx
{
  namespace bcond
  {
    using namespace arakawa_c;

// TODO: add detail namespace, move to detail directory?
    template <typename real_t>
    class cyclic_left_common : public bcond_t<real_t>
    {
      protected:

      // member fields
      rng_t 
        left_halo_sclr, rght_edge_sclr,
        left_halo_vctr, rght_edge_vctr;
      const int halo;

      public:

      // ctor
      cyclic_left_common(const rng_t &i, const int halo) :
        halo(halo),
        // sclr
	left_halo_sclr(i.first() - halo    , i.first() - 1), // TODO: less repetitions!
	rght_edge_sclr(i.last()  - halo + 1, i.last()     ), // TODO: less repetitions!
        // vector
        left_halo_vctr(
          halo == 1
            ? rng_t::all() // there's no vector halo for halo=1
            : rng_t((i-h).first() - (halo-1)    , (i-h).first() - 1)
        ), // TODO: less repetitions!
        rght_edge_vctr(
          halo == 1
            ? rng_t::all() // there's no vector halo for halo=1
            : rng_t((i+h).last()  - (halo-1) + 1, (i+h).last()     )  // TODO: less repetitions!
        )
      {} 
    };

    template <typename real_t>
    class cyclic_rght_common : public bcond_t<real_t>
    {
      protected:

      // member fields
      rng_t 
        left_edge_sclr, rght_halo_sclr,
        left_edge_vctr, rght_halo_vctr;
      const int halo;

      public:

      // ctor
      cyclic_rght_common(const rng_t &i, const int halo) :
        halo(halo),
        // sclr
	rght_halo_sclr(i.last()  + 1       , i.last()  + halo    ), // TODO: less repetitions!
	left_edge_sclr(i.first()           , i.first() + halo - 1), // TODO: less repetitions!
        // vctr
	rght_halo_vctr(
          halo == 1
            ? rng_t::all() // there's no vector halo for halo=1
            : rng_t((i+h).last()  + 1       , (i+h).last()  + (halo-1)    )
        ), // TODO: less repetitions!
	left_edge_vctr(
          halo == 1
            ? rng_t::all() // there's no vector halo for halo=1
            : rng_t((i-h).first()           , (i-h).first() + (halo-1) - 1)  // TODO: less repetitions!
        )
      {} 
    };
  }; // namespace bcond
}; // namespace libmpdataxx
