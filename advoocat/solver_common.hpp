/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <vector>

namespace solvers
{
  template <int n_eqs, typename _real_t>
  class solver_common
  {
    protected: 

    std::vector<int> n;

    // helper methods invoked by solve()
    virtual void advop(int e) = 0;
    void advop_all()
    {
      for (int e = 0; e < n_eqs; ++e) advop(e);
    }

    void cycle(int e) 
    { 
      n[e] = (n[e] + 1) % 2 - 2; 
    }
    void cycle_all()
    { 
      for (int e = 0; e < n_eqs; ++e) cycle(e);
    }

    virtual void xchng(int e) = 0;
    void xchng_all() 
    {   
      for (int e = 0; e < n_eqs; ++e) xchng(e);
    }

    public:

    typedef _real_t real_t;
 
    solver_common() :
      n(n_eqs, 0) 
    {}

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
}; // namespace solvers
