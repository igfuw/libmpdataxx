/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/concurr/detail/concurr_common.hpp>

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
    class serial : public detail::concurr_common<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
    {
      using parent_t = detail::concurr_common<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
 

      struct mem_t : parent_t::mem_t
      {
        int rank() { return 0; }

	static int size() { return 1; }

        void barrier() { }

        // ctors
        mem_t(const std::array<int, solver_t::n_dims> &span) 
          : parent_t::mem_t(span, size()) 
        {};
      };

      void solve(int nt)
      {
        this->algos[0].solve(nt);
      }

      public:

      // ctor
      serial(const typename solver_t::rt_params_t &p) : 
        parent_t(p, new mem_t(p.span), mem_t::size())
      {}

    };
  }; // namespace concurr
}; // namespace libmpdataxx
