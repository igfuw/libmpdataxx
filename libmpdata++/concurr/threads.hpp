/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// TODO: rename to auto?

#pragma once

#ifdef _OPENMP
#  include <libmpdata++/concurr/openmp.hpp>
#else
#  include <libmpdata++/concurr/boost_thread.hpp>
#endif

namespace libmpdataxx
{
  namespace concurr
  {
    /// @brief shared-memory concurency logic using threads 
    ///        (\ref libmpdataxx::concurr::openmp if supported, 
    ///        \ref libmpdataxx::concurr::boost_thread otherwise)
    template <
      class solver_t,
      bcond::bcond_e bcx,
      bcond::bcond_e bcy = bcond::null,
      bcond::bcond_e bcz = bcond::null
    > using threads = 
#if defined(_OPENMP)
    openmp<solver_t, bcx, bcy, bcz>;
#else
    boost_thread<solver_t, bcx, bcy, bcz>;
#endif
  }; // namespace concurr
}; // namespace libmpdataxx
