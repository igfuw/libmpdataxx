// common code for polar boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/bcond/detail/bcond_common.hpp>
#include <libmpdata++/formulae/arakawa_c.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    namespace detail
    {
      using namespace arakawa_c;

      template <typename real_t, int halo, int n_dims>
      class polar_common : public bcond_common<real_t, halo, n_dims>
      {
        using parent_t = bcond_common<real_t, halo, n_dims>;

        protected:

        // member fields
        const int pole;

        int polar_neighbours(const int j)
        {
          return (j + pole) % (2 * pole);
        }

        public:

        // ctor
        polar_common(
          const rng_t &i,
          const std::array<int, n_dims> &grid_size
        ) :
          parent_t(i, grid_size),
          pole((grid_size[0] - 1) / 2)
        {}
      };
    } // namespace detail
  } // namespace bcond
} // namespace libmpdataxx
