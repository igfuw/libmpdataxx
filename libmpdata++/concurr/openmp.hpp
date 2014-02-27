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
        mem_t(const std::array<int, solver_t::n_dims> &span) : parent_t::mem_t(span, size()) {};
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

      // ctor
      openmp(const typename solver_t::rt_params_t &p) : 
        parent_t(p, new mem_t(p.span), mem_t::size())
      {}

    };
  }; // namespace concurr
}; // namespace libmpdataxx
