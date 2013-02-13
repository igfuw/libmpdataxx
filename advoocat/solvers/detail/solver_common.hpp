/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "../../blitz.hpp"
#include "../../concurr/detail/sharedmem.hpp" 

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      constexpr int max(const int a, const int b)
      {
        return a > b ? a : b;
      }

      template <typename real_t_, int n_dims_, int n_eqs_, int n_tlev, int halo_>
      class solver_common
      {
	public:

        // using enums as "static const int" would need instantiation
        enum { halo = halo_ }; 
        enum { n_dims = n_dims_ };
        enum { n_eqs = n_eqs_ };

        typedef real_t_ real_t;
        typedef blitz::Array<real_t_, n_dims_> arr_t;

	protected: 

        typedef concurr::detail::sharedmem<real_t, n_dims, n_eqs> mem_t;
	mem_t *mem;

	// helper methods invoked by solve()
	virtual void advop(int e) = 0;
	void advop_all()
	{
	  for (int e = 0; e < n_eqs; ++e) advop(e);
	}

	void cycle(int e) 
	{ 
	  mem->n[e] = (mem->n[e] + 1) % n_tlev - n_tlev;  // TODO: - n_tlev not needed?
	}
	void cycle_all()
	{ 
	  for (int e = 0; e < n_eqs; ++e) cycle(e);
	}

	virtual void xchng(int e, int l = 0) = 0;
	void xchng_all() 
	{   
	  for (int e = 0; e < n_eqs; ++e) xchng(e);
	}

	public:
      
	// ctor
	solver_common(mem_t *mem) :
	  mem(mem)
	{ }

	void solve(const int nt) 
	{   
	  for (int t = 0; t < nt; ++t) 
	  {   
	    xchng_all();
	    advop_all();
	    cycle_all();
	  }   
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
