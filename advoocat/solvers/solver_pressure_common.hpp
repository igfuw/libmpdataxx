/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "../formulae/courant_formulae.hpp"

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
	typedef typename parent_t::mem_t mem_t;
	typedef typename parent_t::real_t real_t;

	struct params_t : parent_t::params_t { };

	// member fields
	rng_t im, jm;
	real_t dx = 1, dz = 1;  //TODO don't assume dx=dz=1

        virtual void ini_solver() = 0;

	void ini_courant()
	{
	  this->xchng(u);
	  this->xchng(w);

	  formulae::courant::intrp<0>(this->mem.C[0], this->psi(u), im, this->j^this->halo, this->dt, dx);
	  formulae::courant::intrp<1>(this->mem.C[1], this->psi(u), jm, this->i^this->halo, this->dt, dz);
	}

	void update_courant()
	{
	  extrp_velocity(u);      //extrapolate velocity field in time (t+1/2)
	  extrp_velocity(w);

	  this->xchng(u, 1);      // filling halos for velocity filed
	  this->xchng(w, 1);      // psi[n-1] was overwriten for that by extrp_velocity

	  formulae::courant::intrp<0>(this->mem.C[0], this->psi(u, -1), im, this->j^this->halo, this->dt, dx);
	  formulae::courant::intrp<1>(this->mem.C[1], this->psi(w, -1), jm, this->i^this->halo, this->dt, dz);
	}

	void extrp_velocity(int e) // extrapolate in time to t+1/2
	{            // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
	  auto tmp = this->psi(e, -1);
	  tmp /= -2;
	  tmp += 3./2 * this->psi(e);
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
              this->ini_solver();
	      pressure_solver_update(this->dt);
	      pressure_solver_apply(this->dt);
	    }
	    if (t!=0)
	    {
	      std::cerr<<"t= "<<t<<std::endl;
	      update_courant();
	      this->forcings(this->dt / 2);
	      pressure_solver_apply(this->dt);
	      inhomo_solver_t::parent_t::solve(1);
	      // this->xchng_all();
	      // this->advop_all();
	      // this->cycle_all(); 
	      this->forcings(this->dt / 2);
	      pressure_solver_update(this->dt);
	      pressure_solver_apply(this->dt);
	    }
	  }
	  //std::cerr<<"total number of pseudotime iterations = "<<iters<<std::endl;
	}

	// ctor
	pressure_solver_common(
	  mem_t &mem,
	  const rng_t &i, 
	  const rng_t &j, 
	  const params_t &p
	) : 
	  parent_t(mem, i, j, p),
	  im(i.first() - 1, i.last()),
	  jm(j.first() - 1, j.last())
	{} 

	static void alloc(mem_t &mem, const int nx, const int ny)
	{
	  parent_t::alloc(mem, nx, ny);
          // TODO: for sure there will be some common temporaries!
        }
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
