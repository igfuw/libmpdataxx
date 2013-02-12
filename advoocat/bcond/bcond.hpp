/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

namespace advoocat
{
  namespace bcond
  {
    enum bcond_e { null, cyclic }; 

    template <typename real_t>
    class bcond_t
    {
      public:

      // 1D
      virtual void fill_halos(const blitz::Array<real_t, 1> &) 
      {
        assert(false && "bcond::fill_halos() called!");
      };

      // 2D
      virtual void fill_halos(const blitz::Array<real_t, 2> &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos() called!");
      };

      // 3D
      virtual void fill_halos(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos() called!");
      };
    };
  }; // namespace bcond
}; // namespace advoocat
