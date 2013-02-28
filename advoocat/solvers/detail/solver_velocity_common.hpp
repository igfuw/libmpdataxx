/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "../../formulae/courant_formulae.hpp"

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      template <class inhomo_solver_t, int u, int w>
      class solver_velocity_common : public inhomo_solver_t
      {
	protected:

	using parent_t = inhomo_solver_t;
	typedef typename parent_t::real_t real_t;

	// member fields
        const rng_t im, jm;
	const real_t dx, dz;

	void hook_ante_loop()
	{
          rng_t &i = this->i, &j = this->j;
          // allow extrapolation at the first time-step
          this->state(u, -1)(i, j) = this->state(u)(i, j);
          this->state(w, -1)(i, j) = this->state(w)(i, j);
	}

	void extrp_velocity(int e) // extrapolate in time to t+1/2
	{            // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
          rng_t &i = this->i, &j = this->j;
	  auto tmp = this->state(e, -1);

	  tmp(i,j) /= -2;
	  tmp(i,j) += 3./2 * this->state(e)(i,j);
	}

	void hook_ante_step()
	{
	  extrp_velocity(u);      //extrapolate velocity field in time (t+1/2)
	  extrp_velocity(w);

	  this->xchng(u, -1);      // filling halos for velocity filed
	  this->xchng(w, -1);      // psi[n-1] was overwriten for that by extrp_velocity

	  formulae::courant::intrp<0>(this->mem->C[0], this->state(u, -1), im, this->j^this->halo, this->dt, dx);
	  formulae::courant::intrp<1>(this->mem->C[1], this->state(w, -1), jm, this->i^this->halo, this->dt, dz);

          parent_t::hook_ante_step(); // forcings
	}

	public:

	struct params_t : parent_t::params_t 
        { 
          real_t dx, dz;
        };

	// ctor
	solver_velocity_common(
	  typename parent_t::mem_t *mem,
          typename parent_t::bc_p &bcxl,
          typename parent_t::bc_p &bcxr,
          typename parent_t::bc_p &bcyl,
          typename parent_t::bc_p &bcyr,
	  const rng_t &i, 
	  const rng_t &j, 
	  const params_t &p
	) : 
	  parent_t(mem, bcxl, bcxr, bcyl, bcyr, i, j, p),
          im(i.first() - 1, i.last()),
          jm(j.first() - 1, j.last()),
          dx(p.dx), dz(p.dz)
	{} 
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
