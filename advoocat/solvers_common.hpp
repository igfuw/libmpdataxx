/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

namespace solvers
{
  template <int n_eqs>
  class solver_common
  {
    protected: 

    std::vector<int> n;

    virtual void advop(int e) = 0;

    // helper methods invoked by solve()
    void adv()
    {
      for (int e = 0; e < n_eqs; ++e) advop(e);
    }

    void cycle(int e = -1) 
    { 
      if (e == -1) 
        for (e = 0; e < n_eqs; ++e) cycle(e);
      else
        n[e] = (n[e] + 1) % 2 - 2; 
    }

    public:
 
    solver_common() :
      n(n_eqs, 0) 
    {}
  };
}; // namespace solvers
