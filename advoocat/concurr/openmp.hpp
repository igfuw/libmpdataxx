/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "concurr_common.hpp"

#ifdef _OPENMP
# include <omp.h>
#endif

namespace advoocat
{
  namespace concurr
  {
    template <
      class solver_t,
      bcond::bcond_e bcx,
      bcond::bcond_e bcy = bcond::null,
      bcond::bcond_e bcz = bcond::null
    >
    class openmp : public detail::concurr_common<solver_t, bcx, bcy, bcz>
    {
      using parent_t = detail::concurr_common<solver_t, bcx, bcy, bcz>;
 
      int size() 
      {
#if defined(_OPENMP)
        return omp_get_num_procs();
#else
        return 1;
#endif
      }

      struct mem_t : parent_t::mem_t
      {
        void barrier()
        {
#pragma omp barrier
        }
      };

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
      ) : 
        parent_t(s0, params, size())
      {}

      // 2D ctor
      openmp(
	const int s0,
	const int s1,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) : 
	parent_t(s0, s1, params, size(), 1)
      {}

      // 3D ctor
      openmp(
	const int s0,
	const int s1,
	const int s2,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) :
	parent_t(s0, s1, s2, params, size(), 1, 1)
      {}
    };
  }; // namespace concurr
}; // namespace advoocat
