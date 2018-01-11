/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/concurr/detail/concurr_common.hpp>

#include <limits>

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
	static int size(const unsigned max_threads = std::numeric_limits<unsigned>::max())
	{
#if defined(_OPENMP)
	  const char *env_var("OMP_NUM_THREADS");

	  int nthreads = std::min(max_threads, static_cast<unsigned>(
            (std::getenv(env_var) != NULL) ?  std::atoi(std::getenv(env_var)) : omp_get_max_threads()
          ));

          omp_set_num_threads(nthreads);

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
        mem_t(const std::array<int, solver_t::n_dims> &grid_size) : parent_t::mem_t(grid_size, size(grid_size[0])) {};
      };

      void solve(typename parent_t::advance_arg_t nt)
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
        parent_t(p, new mem_t(p.grid_size), mem_t::size(p.grid_size[0]))
      {}

    };
  } // namespace concurr
} // namespace libmpdataxx
