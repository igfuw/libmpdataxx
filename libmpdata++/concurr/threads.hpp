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
      bcond::bcond_e bcxl,
      bcond::bcond_e bcxr,
      bcond::bcond_e bcyl = bcond::null,
      bcond::bcond_e bcyr = bcond::null,
      bcond::bcond_e bczl = bcond::null,
      bcond::bcond_e bczr = bcond::null
    > using threads = 
#if defined(_OPENMP)
    openmp<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
#else
    boost_thread<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
#endif
  } // namespace concurr
} // namespace libmpdataxx
