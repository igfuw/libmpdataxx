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
      template <class inhomo_solver_t, int u, int w, int tht>
      class pressure_solver_common : public inhomo_solver_t
      {
	protected:

	using parent_t = inhomo_solver_t;
	typedef typename parent_t::real_t real_t;

	struct params_t : parent_t::params_t { };

	// member fields
	rng_t im, jm;
	real_t dx = 1, dz = 1;  //TODO don't assume dx=dz=1
        int iters = 0;

	void ini_courant()
	{
	  this->xchng(u);
	  this->xchng(w);

	  formulae::courant::intrp<0>(this->mem->C[0], this->state(u), im, this->j^this->halo, this->dt, dx);
	  formulae::courant::intrp<1>(this->mem->C[1], this->state(w), jm, this->i^this->halo, this->dt, dz);
	}

	void extrp_courant()
	{
	  extrp_velocity(u);      //extrapolate velocity field in time (t+1/2)
	  extrp_velocity(w);

	  this->xchng(u, -1);      // filling halos for velocity filed
	  this->xchng(w, -1);      // psi[n-1] was overwriten for that by extrp_velocity

	  formulae::courant::intrp<0>(this->mem->C[0], this->state(u, -1), im, this->j^this->halo, this->dt, dx);
	  formulae::courant::intrp<1>(this->mem->C[1], this->state(w, -1), jm, this->i^this->halo, this->dt, dz);
	}

	void extrp_velocity(int e) // extrapolate in time to t+1/2
	{            // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
          rng_t &i = this->i, &j = this->j;
	  auto tmp = this->state(e, -1);

	  tmp(i,j) /= -2;
	  tmp(i,j) += 3./2 * this->state(e)(i,j);
	}

	virtual void ini_pressure() = 0;
	virtual void pressure_solver_update(real_t dt) = 0;
	virtual void pressure_solver_apply(real_t dt) = 0;

	public:

	void solve(int nt) 
	{   
	  for (int t = 0; t < nt; ++t)
	  {   
	    if (t==0)
	    {
	      ini_courant();
	      ini_pressure();
	      this->forcings(this->dt / 2);
	      inhomo_solver_t::parent_t::solve(1);
	      this->forcings(this->dt / 2);
	      pressure_solver_update(this->dt);
// Tu równolegle state(u) jest losowe
std::ostringstream s;
s << "rank=" << this->mem->rank() 
  << " sum(state(u)) = " << this->mem->sum(this->state(u)(this->i, this->j)) << std::endl;
	      pressure_solver_apply(this->dt);
s << " sum(state(u)) = " << this->mem->sum(this->state(u)(this->i, this->j)) << std::endl;
std::cerr << s.str();
	    }
	    if (t!=0)
	    {
	      std::cerr<<"t= "<<t<<std::endl;
	      extrp_courant();
	      this->forcings(this->dt / 2);
	      pressure_solver_apply(this->dt);
	      inhomo_solver_t::parent_t::solve(1);
	      this->forcings(this->dt / 2);
	      pressure_solver_update(this->dt);
	      pressure_solver_apply(this->dt);
	    }
	  }
std::cerr<<"number of pseudo time iterations "<<iters<<std::endl;
	}

	// ctor
	pressure_solver_common(
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
	  jm(j.first() - 1, j.last())
	{} 

	static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
	{
	  parent_t::alloc(mem, nx, ny);
          // TODO: for sure there will be some common temporaries!
        }
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
