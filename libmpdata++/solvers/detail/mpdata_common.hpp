/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

//#include <array>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template<typename real_t, int n_iters_>
      class mpdata_common 
      {
	protected:

	static const int n_iters = n_iters_;
	static_assert(n_iters > 0, "n_iters <= 0");
      };
    }; // namespace detail
  }; // namespace solvers
}; // namescpae libmpdataxx
