/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

namespace libmpdataxx
{
  namespace formulae
  {
    template <class a_t>
    inline void kahan_zro(a_t c, const a_t&, const a_t&, a_t sum)
    {
      sum = 0;
      c = 0;
    }

    template <class a_t, class f_t>
    inline void kahan_add(a_t c, a_t y, a_t t, a_t sum, f_t input)
    {
      y = input - c;
      t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
  }
}
