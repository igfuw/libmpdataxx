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

      rng_t west_polar_halo(const rng_t &i)
      {
        return rng_t(i.first(), i.first() + this->halo + pole);
      }

      rng_t east_polar_halo(const rng_t &i)
      {
        return rng_t(i.first() + this->halo + pole + 1, i.last());
      }

      rng_t west_sclr_polar_edge(const rng_t &i)
      {
        return rng_t(i.first() + this->halo + 1, i.first() + 2 * this->halo + pole);
      }

      rng_t east_sclr_polar_edge(const rng_t &i)
      {
        return rng_t(i.first() + pole, i.last() - this->halo);
      }
      
      rng_t west_vctr_polar_edge(const rng_t &i)
      {
        return rng_t(i.first() + this->halo + 1, i.first() + 2 * this->halo + pole - 1);
      }

      rng_t east_vctr_polar_edge(const rng_t &i)
      {
        return rng_t(i.first() + pole + 1, i.last() - this->halo);
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
