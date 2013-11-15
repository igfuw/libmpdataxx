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
      bcond::bcond_e bcx,
      bcond::bcond_e bcy = bcond::null,
      bcond::bcond_e bcz = bcond::null
    >
    class serial : public detail::concurr_common<solver_t, bcx, bcy, bcz>
    {
      using parent_t = detail::concurr_common<solver_t, bcx, bcy, bcz>;
 

      struct mem_t : parent_t::mem_t
      {
        int rank() { return 0; }

	static int size() { return 1; }

        void barrier() { }

        // ctors
        mem_t(const int n_eqs, const int s0) : parent_t::mem_t(n_eqs, s0, size()) {};
        mem_t(const int n_eqs, const int s0, const int s1) : parent_t::mem_t(n_eqs, s0, s1, size()) {};
        mem_t(const int n_eqs, const int s0, const int s1, const int s2) : parent_t::mem_t(n_eqs, s0, s1, s2, size()) {};
      };

      void solve(int nt)
      {
        this->algos[0].solve(nt);
      }

      public:

// TODO: this is an exact duplicate from openmp :(
      // 1D ctor
      serial(
	const int s0,
	const typename solver_t::params_t &p = typename solver_t::params_t()
      ) : 
        parent_t(s0, p, new mem_t(p.n_eqs, s0), mem_t::size())
      {}

      // 2D ctor
      serial(
	const int s0,
	const int s1,
	const typename solver_t::params_t &p = typename solver_t::params_t()
      ) : 
	parent_t(s0, s1, p, new mem_t(p.n_eqs, s0, s1), mem_t::size(), 1)
      {}

      // 3D ctor
      serial(
	const int s0,
	const int s1,
	const int s2,
	const typename solver_t::params_t &p = typename solver_t::params_t()
      ) :
	parent_t(s0, s1, s2, p, new mem_t(p.n_eqs, s0, s1, s2), mem_t::size(), 1, 1)
      {}
    };
  }; // namespace concurr
}; // namespace libmpdataxx
