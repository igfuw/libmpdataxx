/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "../formulae/mpdata_formulae.hpp"
#include "../formulae/donorcell_formulae.hpp"
#include "solver_2d.hpp"

namespace advoocat
{
  namespace solvers
  {
    using namespace arakawa_c;

    template<int n_iters, class bcx_t, class bcy_t, class mem_t>
    class mpdata_2d : public detail::solver_2d<bcx_t, bcy_t, mem_t>
    {
      static_assert(n_iters > 0, "n_iters <= 0");

      using parent_t = detail::solver_2d<bcx_t, bcy_t, mem_t>;
      using arr_2d_t = typename mem_t::arr_t;

      static const int n_tmp = n_iters > 2 ? 2 : 1;

      // member fields
      arrvec_t<arr_2d_t> *tmp[2];
      rng_t im, jm;

      protected:

      // method invoked by the solver
      void advop(int e)
      {
	for (int step = 0; step < n_iters; ++step) 
	{
	  if (step == 0) 
	    formulae::donorcell::op_2d(
	      this->mem.psi[e], 
	      this->mem.n[e], this->mem.C, this->i, this->j);
	  else
	  {
	    this->cycle(e);
	    this->bcx.fill_halos(this->mem.psi[e][this->mem.n[e]], this->j^this->halo);
	    this->bcy.fill_halos(this->mem.psi[e][this->mem.n[e]], this->i^this->halo);

	    // choosing input/output for antidiff C
            const arrvec_t<arr_2d_t>
	      &C_unco = (step == 1) 
		? this->mem.C 
		: (step % 2) 
		  ? *tmp[1]  // odd steps
		  : *tmp[0], // even steps
	      &C_corr = (step  % 2) 
		? *tmp[0]    // odd steps
		: *tmp[1];   // even steps

	    // calculating the antidiffusive C 
	    C_corr[0](this->im+h, this->j) = 
	      formulae::mpdata::antidiff<0>(
		this->mem.psi[e][this->mem.n[e]], 
		this->im, this->j, C_unco
	      );
	    this->bcy.fill_halos(C_corr[0], this->i^h);

	    C_corr[1](this->i, this->jm+h) = 
	      formulae::mpdata::antidiff<1>(
              this->mem.psi[e][this->mem.n[e]], 
              this->jm, this->i, C_unco
	    );
	    this->bcx.fill_halos(C_corr[1], this->j^h);

	    // donor-cell step 
	    formulae::donorcell::op_2d(this->mem.psi[e], 
	      this->mem.n[e], C_corr, this->i, this->j);
	  }
	}
      }

      public:

      struct params_t {};

      // ctor
      mpdata_2d(mem_t &mem, const rng_t &i, const rng_t &j, const params_t &) : 
	detail::solver_2d<bcx_t, bcy_t, mem_t>(mem, i, j, 1), 
	im(i.first() - 1, i.last()),             //TODO get correct version from sylwester
	jm(j.first() - 1, j.last())
      {
	for (int n = 0; n < n_tmp; ++n)
          tmp[n] = &mem.tmp[std::string(__FILE__)][n];
      }

      // memory allocation (to be called once per shared-mem node)
      static void alloctmp(
        std::unordered_map<std::string, boost::ptr_vector<arrvec_t<arr_2d_t>>> &tmp, 
        const int nx,
        const int ny
      )   
      {   
        const rng_t i(0, nx-1), j(0, ny-1);
        const int halo = 1;  //TODO get correct version from sylwester
        for (int n = 0; n < n_tmp; ++n)
        {
          tmp[std::string(__FILE__)].push_back(new arrvec_t<arr_2d_t>());
          tmp[std::string(__FILE__)].back().push_back(new arr_2d_t( i^h, j^halo ));
          tmp[std::string(__FILE__)].back().push_back(new arr_2d_t( i^halo, j^h ));
        }
      }   
    };
  }; // namespace solvers
}; // namespace advoocat
