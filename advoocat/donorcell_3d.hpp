/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_3d.hpp"
#include "donorcell_formulae.hpp"

namespace solvers
{
  template<
    class bcx_t, 
    class bcy_t, 
    class bcz_t, 
    class mem_t
  > 
  class donorcell_3d : public solver_3d<bcx_t, bcy_t, bcz_t, mem_t> 
  {
    void advop(int e)
    {
      donorcell::op_3d(
        this->mem.psi[e], this->mem.n[e], this->mem.C, this->i, this->j, this->k
      );
    }

    public:

    struct params {};

    // ctor
    donorcell_3d(
      mem_t &mem, 
      const rng_t &i, 
      const rng_t &j, 
      const rng_t &k, 
      const params &
    ) :
      solver_3d<bcx_t, bcy_t, bcz_t, mem_t>(mem, i, j, k, /* halo = */ 1)
    {}  
  };
}; // namespace solvers

