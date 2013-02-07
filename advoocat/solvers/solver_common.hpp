/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "../sharedmem.hpp"

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      template <class sharedmem>
      class solver_common
      {
	public:

	typedef sharedmem mem_t;

	protected: 

	mem_t &mem;

	const int n_tlev = 2;

	// helper methods invoked by solve()
	virtual void advop(int e) = 0;
	void advop_all()
	{
	  for (int e = 0; e < mem_t::n_eqs; ++e) advop(e);
	}

	void cycle(int e) 
	{ 
	  mem.n[e] = (mem.n[e] + 1) % n_tlev - n_tlev; 
	}
	void cycle_all()
	{ 
	  for (int e = 0; e < mem_t::n_eqs; ++e) cycle(e);
	}

	virtual void xchng(int e) = 0;
	void xchng_all() 
	{   
	  for (int e = 0; e < mem_t::n_eqs; ++e) xchng(e);
	}

	public:
      
	// ctor
	solver_common(mem_t &mem) :
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
