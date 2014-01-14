/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/solvers/detail/solver_1d.hpp>
#include <libmpdata++/formulae/donorcell_formulae.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    template<
      typename ct_params_t, 
      int halo = formulae::donorcell::halo
    >
    class donorcell_1d : public detail::solver<
      ct_params_t,
      1,
      formulae::donorcell::n_tlev, 
      detail::max(halo, formulae::donorcell::halo)
    > 
    {
      using parent_t = detail::solver<
        ct_params_t,
        1,
        formulae::donorcell::n_tlev, 
        detail::max(halo, formulae::donorcell::halo)
      >;
   
      void advop(int e)
      {
        formulae::donorcell::op_1d<ct_params_t::opts>(
          this->mem->psi[e], 
	  this->mem->GC[0], 
	  *this->mem->G, 
	  this->n[e], 
	  this->i
        );
      }

      public:

      // ctor
      donorcell_1d(
        typename parent_t::ctor_args_t args, 
        const typename parent_t::params_t &p
      ) :
        parent_t(args, p)
      {}  
    };
  }; // namespace solvers
}; // namespace libmpdataxx
