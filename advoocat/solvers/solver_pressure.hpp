/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include "solver_common.hpp"
#include "solver_inhomo.hpp"
#include "../formulae/courant_formulae.hpp"
#include "../formulae/nabla_formulae.hpp" //gradient, diveregnce

namespace advoocat
{
  namespace solvers
  {
    template <class parent_t_, int u, int w, int tht>
    class pressure_solver : public parent_t_
    {
      public:

      using parent_t = parent_t_;
      typedef typename parent_t::mem_t mem_t;
      typedef typename parent_t::real_t real_t;
      using arr_2d_t = typename mem_t::arr_t;

      real_t Prs_amb;

      //TODO probably don't need those
      arr_2d_t Phi, tmp_u, tmp_w, tmp_div, tmp_x, tmp_z;

      virtual void forcings(real_t dt) = 0;

      private:

      rng_t im, jm;
      //TODO don't assume dx=dz=1
      real_t dx = 1;
      real_t dz = 1;

      void extrp_velocity(int e) // extrapolate in time to t+1/2
      {            // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
	auto tmp = this->psi(e, -1);
	tmp /= -2;
	tmp += 3./2 * this->psi(e);
      }

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

      void ini_pressure()
      {
	// dt/2 * (Prs-Prs_amb) / rho
	Phi(this->i^this->halo, this->j^this->halo) = real_t(0);
      }

      void pressure_solver_update(real_t dt)
      {
	using namespace arakawa_c;
	using formulae::nabla_op::grad;
	using formulae::nabla_op::div;

	real_t beta = .25;  //TODO
	real_t rho = 1.;   //TODO    

	int halo = this->halo;
	rng_t i = this->i;
	rng_t j = this->j;

	tmp_u = this->psi(u);
	tmp_w = this->psi(w);

    std::cerr<<"--------------------------------------------------------------"<<std::endl;
	//pseudo-time loop
	real_t err = 1.;
	while (err > .01)
	{
	  this->xchng(Phi,   i^halo, j^halo);
	  this->xchng(tmp_u, i^halo, j^halo);
	  this->xchng(tmp_w, i^halo, j^halo);

	  tmp_x(i, j) = rho * tmp_u(i, j) - grad<0>(Phi(i^halo, j^halo), i, j, real_t(1));
	  tmp_z(i, j) = rho * tmp_w(i, j) - grad<1>(Phi(i^halo, j^halo), j, i, real_t(1));
     
	  this->xchng(tmp_x, i^halo, j^halo);
	  this->xchng(tmp_z, i^halo, j^halo);

	  tmp_div(i, j) = div(tmp_x(i^halo,j^halo), tmp_z(i^halo, j^halo), i, j, real_t(1), real_t(1));
	  Phi(i, j) -= beta * tmp_div(i, j);

	  this->xchng(tmp_u, i^halo, j^halo);
	  this->xchng(tmp_w, i^halo, j^halo);

	  err = std::abs(max(tmp_div)) + std::abs(min(tmp_div));
std::cerr<<"div:  ( "<<min(tmp_div)<<" --> "<<max(tmp_div)<<" )"<<std::endl;
	}
	//end of pseudo_time loop
	this->xchng(this->Phi, i^halo, j^halo);

	tmp_u(i, j) -= grad<0>(Phi(i^halo, j^halo), i, j, real_t(1));
	tmp_w(i, j) -= grad<1>(Phi(i^halo, j^halo), j, i, real_t(1));

	tmp_u -= this->psi(u);
	tmp_w -= this->psi(w);
      }

      void pressure_solver_apply(real_t dt)
      {
	auto U = this->psi(u);
	auto W = this->psi(w);

	U += tmp_u;
	W += tmp_w;
      }

      public:

      struct params_t : parent_t::params_t { real_t Prs_amb; };

      // ctor

      // ctor
      pressure_solver(
	mem_t &mem,
	const rng_t &i,
	const rng_t &j,
	const params_t &p
      ) :
	parent_t(mem, i, j, p),
	Prs_amb(p.Prs_amb),
	tmp_div(mem.tmp[std::string(__FILE__)][0][0]),
	tmp_x(mem.tmp[std::string(__FILE__)][0][1]),
	tmp_z(mem.tmp[std::string(__FILE__)][0][2]),
	tmp_u(mem.tmp[std::string(__FILE__)][0][3]),
	tmp_w(mem.tmp[std::string(__FILE__)][0][4]),
	Phi(mem.tmp[std::string(__FILE__)][0][5]),
	im(i.first() - 1, i.last()),
	jm(j.first() - 1, j.last())
      {}

      static void alloctmp(
        std::unordered_map<std::string, boost::ptr_vector<arrvec_t<arr_2d_t>>> &tmp,  
        const int nx, 
        const int ny
      )
      {
        parent_t::alloctmp(tmp, nx, ny);
        const rng_t i(0, nx-1), j(0, ny-1);
        const int hlo = formulae::mpdata::halo; // TODO!!!
        tmp[std::string(__FILE__)].push_back(new arrvec_t<arr_2d_t>());
        tmp[std::string(__FILE__)].back().push_back(new arr_2d_t( i, j )); 
        for (int n=0; n < 5; ++n) 
          tmp[std::string(__FILE__)].back().push_back(new arr_2d_t( i^hlo, j^hlo )); 
      }

      void solve(int nt)
      {
	for (int t = 0; t < nt; ++t)
	{
	  if (t==0)
	  {
	    ini_courant();
	    ini_pressure();
	    forcings(this->dt / 2);
	    parent_t::parent_t::solve(1);
	    forcings(this->dt / 2);
	    pressure_solver_update(this->dt);
	    pressure_solver_apply(this->dt);
	  }
	  if (t!=0)
	  {
    std::cerr<<"t= "<<t<<std::endl;
	    update_courant();
	    forcings(this->dt / 2);
	    pressure_solver_apply(this->dt);
	    parent_t::parent_t::solve(1);
	    // this->xchng_all();
	    // this->advop_all();
	    // this->cycle_all(); 
	    forcings(this->dt / 2);
	    pressure_solver_update(this->dt);
	    pressure_solver_apply(this->dt);
	  }
	}
      }
    }; 
  }; // namespace solvers
}; // namespace advoocat
