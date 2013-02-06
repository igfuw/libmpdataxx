/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "equip.hpp"

#ifdef _OPENMP
# include <omp.h>
#endif

namespace advoocat
{
  template <class solver_t>
  class openmp : public equip<solver_t>
  {
    void barrier()
    {
#pragma omp barrier
    }

    public:

    void advance(int nt)
    {
      int i = 0;
#pragma omp parallel private(i)
      {
#if defined(_OPENMP)
        i = omp_get_thread_num();
#endif
        this->algos[i].solve(nt);
      } 
    }

    // 1D ctor
    openmp(
      const int s0,
      const typename solver_t::params_t &params = typename solver_t::params_t()
    )
      : equip<solver_t>(s0, params)
    {}

    // 2D ctor
    openmp(
      const int s0,
      const int s1,
      const typename solver_t::params_t &params = typename solver_t::params_t()
    )
      : equip<solver_t>(s0, s1, params)
    {}

    // 3D ctor
    openmp(
      const int s0,
      const int s1,
      const int s2,
      const typename solver_t::params_t &params = typename solver_t::params_t()
    )
      : equip<solver_t>(s0, s1, s2, params)
    {}
  };
}; // namespace advoocat
