/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/concurr/detail/concurr_common.hpp>

#ifdef _OPENMP
# include <omp.h>
#endif

namespace libmpdataxx
{
  namespace concurr
  {
    template <
      class solver_t,
      bcond::bcond_e bcxl,
      bcond::bcond_e bcxr,
      bcond::bcond_e bcyl = bcond::null,
      bcond::bcond_e bcyr = bcond::null,
      bcond::bcond_e bczl = bcond::null,
      bcond::bcond_e bczr = bcond::null
    >
    class openmp : public detail::concurr_common<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
    {
      using parent_t = detail::concurr_common<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
 

      struct mem_t : parent_t::mem_t
      {
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
          // TODO: if (size() != 1) ???
#pragma omp barrier
        }

        // ctors
        mem_t(const std::array<int, solver_t::n_dims> &grid_size) : parent_t::mem_t(grid_size, size()) {};
      };

      void solve(typename parent_t::real_t tshift)
      {
        int i = 0;
#pragma omp parallel private(i)
        {
#if defined(_OPENMP)
          i = omp_get_thread_num();
#endif
          this->algos[i].solve(tshift);
        } 
      }

      public:

      // ctor
      openmp(const typename solver_t::rt_params_t &p) : 
        parent_t(p, new mem_t(p.grid_size), mem_t::size())
      {}

    };
  }; // namespace concurr
}; // namespace libmpdataxx
