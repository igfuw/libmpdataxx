/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "detail/concurr_common.hpp"

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
 

      struct mem_t : parent_t::mem_t
      {
        int rank()
        {
#if defined(_OPENMP)
          return omp_get_thread_num();
#else
          return 0;
#endif
        }

	static int size() 
	{
#if defined(_OPENMP)
	  return omp_get_max_threads();
#else
	  return 1;
#endif
	}

        void barrier()
        {
#pragma omp barrier
        }

        // ctors
        mem_t(int s0) : parent_t::mem_t(s0, size()) {};
        mem_t(int s0, int s1) : parent_t::mem_t(s0, s1, size()) {};
        mem_t(int s0, int s1, int s2) : parent_t::mem_t(s0, s1, s2, size()) {};
      };

      void solve(int nt)
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

      public:

// TODO: coud it be just one ctor with int[solver_t::n_dims]?
// TODO: document that current;y paralllisation only in one dimension
      // 1D ctor
      openmp(
	const int s0,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) : 
        parent_t(s0, params, new mem_t(s0), mem_t::size())
      {}

      // 2D ctor
      openmp(
	const int s0,
	const int s1,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) : 
	parent_t(s0, s1, params, new mem_t(s0, s1), mem_t::size(), 1)
      {}

      // 3D ctor
      openmp(
	const int s0,
	const int s1,
	const int s2,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) :
	parent_t(s0, s1, s2, params, new mem_t(s0, s1, s2), mem_t::size(), 1, 1)
      {}
    };
  }; // namespace concurr
}; // namespace advoocat
