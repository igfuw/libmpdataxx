/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

namespace advoocat
{
  template <class solver_t>
  class equip
  {
    const int n0 = 1, n1 = 1, n2 = 1;

    int min(int span, int rank, int size) 
    { 
      return rank * span / size; 
    }

    int max(int span, int rank, int size) 
    { 
      return min(span, rank + 1, size) - 1; 
    }

    // member fields
    typename solver_t::mem_t mem;
   
    protected:

    boost::ptr_vector<solver_t> algos; 

    public:

    typedef typename solver_t::mem_t::real_t real_t; 

    // 1D ctor
    equip(
      const int s0, 
      const typename solver_t::params_t &params 
    )
      : mem(s0)
    {
      solver_t::alloctmp(mem.tmp, s0);
      for (int i0 = 0; i0 < n0; ++i0) 
	algos.push_back(new solver_t(
	  mem, 
	  rng_t( min(s0, i0, n0), max(s0, i0, n0)),
	  params
	));
    }

    // 2D ctor
    equip(
      const int s0, 
      const int s1, 
      const typename solver_t::params_t &params = typename solver_t::params_t()
    )
      : mem(s0, s1)
    {
      solver_t::alloctmp(mem.tmp, s0, s1);
      for (int i0 = 0; i0 < n0; ++i0) 
	for (int i1 = 0; i1 < n1; ++i1) 
	  algos.push_back(new solver_t(
          mem, 
          rng_t( min(s0, i0, n0), max(s0, i0, n0)),
          rng_t( min(s1, i1, n1), max(s1, i1, n1)),
          params
        ));
    }

    // 3D ctor
    equip(
      const int s0, 
      const int s1, 
      const int s2, 
      const typename solver_t::params_t &params = typename solver_t::params_t()
    )
      : mem(s0, s1, s2)
    {
      solver_t::alloctmp(mem.tmp, s0, s1, s2);
      for (int i0 = 0; i0 < n0; ++i0) 
	for (int i1 = 0; i1 < n1; ++i1) 
	  for (int i2 = 0; i2 < n2; ++i2) 
	    algos.push_back(new solver_t(
	      mem, 
	      rng_t( min(s0, i0, n0), max(s0, i0, n0)),
	      rng_t( min(s1, i1, n1), max(s1, i1, n1)),
	      rng_t( min(s2, i2, n2), max(s2, i2, n2)),
	      params
	    ));
    }

    typename solver_t::mem_t::arr_t state(int e = 0)
    {
      return mem.state(e);
    }

    typename solver_t::mem_t::arr_t courant(int d = 0)
    {
      return mem.courant(d);
    }
  };
}; // namespace advoocat
