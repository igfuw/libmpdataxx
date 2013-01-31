/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_1d.hpp"
#include "donorcell_formulae.hpp"

namespace solvers
{
  template<class bcx_t, class mem_t>
  class donorcell_1d : public solver_1d<bcx_t, mem_t> 
  {
    void advop(int e)
    {
      donorcell::op_1d(
        this->mem.psi[e], this->mem.n[e], this->mem.C[0], this->i
      );
    }

    public:

    struct params {};

    // ctor
    donorcell_1d(mem_t &mem, const rng_t &i, const params &) :
      solver_1d<bcx_t, mem_t>(mem, i, /* halo = */ 1)
    { }  
  };
}; // namespace solvers

