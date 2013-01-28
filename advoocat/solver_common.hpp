/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "blitz.hpp"
#include <vector>
#include <blitz/array.h>

namespace solvers
{
  template <int n_eqs, int n_dims, typename _real_t>
  class solver_common
  {
    public:

    typedef _real_t real_t;
    typedef blitz::Array<real_t, n_dims> arr_t;

    protected: 

    const int n_tlev = 2;

    // member fields
    arrvec_t<arr_t> C;

    // psi contains model state including halo region
    arrvec_t<arr_t> psi[n_eqs];

    std::vector<int> n;

    // helper methods invoked by solve()
    virtual void advop(int e) = 0;
    void advop_all()
    {
      for (int e = 0; e < n_eqs; ++e) advop(e);
    }

    void cycle(int e) 
    { 
      n[e] = (n[e] + 1) % n_tlev - n_tlev; 
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

    /**
     * @brief fills all arrays with NaNs if in debug mode
     */
    solver_common() :
      n(n_eqs, 0) 
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

    public:

    // accessor method for the Courant number field
    arr_t courant(int d = 0) 
    {   
      return C[d]; 
    }   
  };
}; // namespace solvers
