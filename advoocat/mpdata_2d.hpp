/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "mpdata_formulae.hpp"
#include "donorcell_formulae.hpp"
#include "solver_2d.hpp"

namespace solvers
{
  template<int n_iters, class bcx_t, class bcy_t, int n_eqs = 1, typename real_t = float>
  class mpdata_2d : public solver_2d<bcx_t, bcy_t, n_eqs, real_t>
  {
    using parent = solver_2d<bcx_t, bcy_t, n_eqs, real_t>;
    using arr_2d_t = typename parent::arr_t;

    // member fields
    arrvec_t<arr_2d_t> tmp[2];
    rng_t im, jm;

    protected:

    // method invoked by the solver
    void advop(int e)
    {
      for (int step = 0; step < n_iters; ++step) 
      {
        if (step == 0) 
          donorcell::op_2d(
            this->psi[e], 
            this->n[e], this->C, this->i, this->j);
        else
        {
          this->cycle(e);
          this->bcx.fill_halos(this->psi[e][this->n[e]], this->j^this->halo);
          this->bcy.fill_halos(this->psi[e][this->n[e]], this->i^this->halo);

          // choosing input/output for antidiff C
          const arrvec_t<arr_2d_t>
            &C_unco = (step == 1) 
              ? this->C 
              : (step % 2) 
                ? tmp[1]  // odd steps
                : tmp[0], // even steps
            &C_corr = (step  % 2) 
              ? tmp[0]    // odd steps
              : tmp[1];   // even steps

          // calculating the antidiffusive C 
          C_corr[0](this->im+h, this->j) = 
            mpdata::antidiff<0>(
              this->psi[e][this->n[e]], 
              this->im, this->j, C_unco
            );
          this->bcy.fill_halos(C_corr[0], this->i^h);

          C_corr[1](this->i, this->jm+h) = 
            mpdata::antidiff<1>(
              this->psi[e][this->n[e]], 
              this->jm, this->i, C_unco
          );
          this->bcx.fill_halos(C_corr[1], this->j^h);

          // donor-cell step 
          donorcell::op_2d(this->psi[e], 
            this->n[e], C_corr, this->i, this->j);
        }
      }
    }

    public:

    // ctor
    mpdata_2d(int nx, int ny) : 
      solver_2d<bcx_t, bcy_t, n_eqs, real_t>(nx, ny, 1), 
      im(this->i.first() - 1, this->i.last()),
      jm(this->j.first() - 1, this->j.last())
    {
      int n_tmp = n_iters > 2 ? 2 : 1;
      for (int n = 0; n < n_tmp; ++n)
      {
        tmp[n].push_back(new arr_2d_t( this->i^h, this->j^this->halo ));
        tmp[n].push_back(new arr_2d_t( this->i^this->halo, this->j^h ));
      }
    }
  };
}; // namespace solvers
