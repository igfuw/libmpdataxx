/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "mpdata_formulae.hpp"
#include "donorcell_formulae.hpp"
#include "solver_1d.hpp"

namespace solvers
{
  template<int n_iters, class bcx_t, class mem_t>
  class mpdata_1d : public solver_1d<bcx_t, mem_t>
  {
    static_assert(n_iters > 0, "n_iters <= 0");

    using parent = solver_1d<bcx_t, mem_t>;
    using arr_1d_t = typename mem_t::arr_t;

    // member fields
    arrvec_t<arr_1d_t> tmp[2];
    rng_t im;

    protected:

    // method invoked by the solver
    void advop(int e)
    {
      for (int step = 0; step < n_iters; ++step) 
      {
        if (step == 0) 
          donorcell::op_1d(this->mem.psi[e], this->mem.n[e], this->mem.C[0], this->i);
        else
        {
          this->cycle(e);
          this->bcx.fill_halos(this->mem.psi[e][this->mem.n[e]]);

          // choosing input/output for antidiff C
          const arr_1d_t
            &C_unco = (step == 1) 
              ? this->mem.C[0] 
              : (step % 2) 
                ? tmp[1][0]  // odd steps
                : tmp[0][0], // even steps
            &C_corr = (step  % 2) 
              ? tmp[0][0]    // odd steps
              : tmp[1][0];   // even steps

          // calculating the antidiffusive C 
          C_corr(im+h) = 
            mpdata::antidiff(
              this->mem.psi[e][this->mem.n[e]], 
              im, C_unco[0]
            );

          // donor-cell step 
          donorcell::op_1d(this->mem.psi[e], 
            this->mem.n[e], C_corr[0], this->i);
        }
      }
    }

    public:

    struct params {};

    // ctor
    mpdata_1d(mem_t &mem, const rng_t &i, const params &) : 
      parent(mem, i, /* halo = */1), 
      im(i.first() - 1, i.last())
    {
      int n_tmp = n_iters > 2 ? 2 : 1;
      for (int n = 0; n < n_tmp; ++n)
        tmp[n].push_back(new arr_1d_t( i^h )); // TODO: alloc()
    }

  };
};
