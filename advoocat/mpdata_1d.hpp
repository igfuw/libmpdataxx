/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "mpdata_common.hpp"
#include "solver_1d.hpp"
#include "donorcell.hpp"

namespace solvers
{
  template<int n_iters, class bcx_t, int n_eqs = 1, typename real_t = float>
  class mpdata_1d : public solver_1d<bcx_t, n_eqs, real_t>
  {
    // member fields
    arrvec_t<arr_1d_t<real_t>> tmp[2];
    rng_t im;

    protected:

    // method invoked by the solver
    void advop(int e)
    {
      for (int step = 0; step < n_iters; ++step) 
      {
        if (step == 0) 
          donorcell::op_1d(this->psi[e], this->n[e], this->C[0], this->i);
        else
        {
          this->cycle(e);
          this->bcx.fill_halos(this->psi[e][this->n[e]]);

          // choosing input/output for antidiff C
          const arr_1d_t<real_t>
            &C_unco = (step == 1) 
              ? this->C[0] 
              : (step % 2) 
                ? tmp[1][0]  // odd steps
                : tmp[0][0], // even steps
            &C_corr = (step  % 2) 
              ? tmp[0][0]    // odd steps
              : tmp[1][0];   // even steps

          // calculating the antidiffusive C 
          C_corr(this->im+h) = 
            mpdata::antidiff(
              this->psi[e][this->n[e]], 
              this->im, C_unco[0]
            );

          // donor-cell step 
          donorcell::op_1d(this->psi[e], 
            this->n[e], C_corr[0], this->i);
        }
      }
    }

    public:

    // ctor
    mpdata_1d(int nx) : 
      solver_1d<bcx_t, n_eqs, real_t>(nx, /* halo = */1), 
      im(this->i.first() - 1, this->i.last())
    {
      int n_tmp = n_iters > 2 ? 2 : 1;
      for (int n = 0; n < n_tmp; ++n)
      {
        tmp[n].push_back(new arr_1d_t<real_t>(
          this->i^h));
      }
    }

  };
};
