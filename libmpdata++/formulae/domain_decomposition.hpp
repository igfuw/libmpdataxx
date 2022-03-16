/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/blitz.hpp>

namespace libmpdataxx
{
  namespace domain_decomposition
  {
    namespace detail
    {
      // helper methods to define subdomain ranges
      int min(const int &span, const int &rank, const int &size)
      {
        return rank * span / size;
      }

      int max(const int &span, const int &rank, const int &size)
      {
        return min(span, rank + 1, size) - 1;
      }
    };

    // get part of 'span' assigned to 'rank' (out of 'size' ranks)
    rng_t slab(
      const rng_t &span,
      const int &rank = 0,
      const int &size = 1
    ) {
      return rng_t(
        span.first() + detail::min(span.length(), rank, size),
        span.first() + detail::max(span.length(), rank, size)
      );
    }
  }
} 
