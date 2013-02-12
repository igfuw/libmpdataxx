/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

/*
eqs 12a and 12b from Smolarkiewicz & Margolin 1994

for it to work one would need another scheme to initialize first guess
for phi[n] and err[n]

also, it's more expansive than the other method
//TODO - finish later(?) or remove
*/


#pragma once
#include "detail/solver_pressure_common.hpp"
#include "../formulae/nabla_formulae.hpp" //gradient, diveregnce

namespace advoocat
{
  namespace solvers
  {
    template <class inhomo_solver_t, int u, int w, int tht>
    class pressure_crk2 : public detail::pressure_solver_common<inhomo_solver_t, u, w, tht>
    {
      public:

      using parent_t = detail::pressure_solver_common<inhomo_solver_t, u, w, tht>;
      typedef typename parent_t::mem_t mem_t;
      typedef typename parent_t::real_t real_t;
      using arr_2d_t = typename mem_t::arr_t;

      int iters = 0;
      
      // timelevels
      int n;

      //TODO probably don't need those
      arr_2d_t tmp_u, tmp_w, tmp_x, tmp_z;
      arr_2d_t lap_err, tmp_e1, tmp_e2;

      arrvec_t<arr_2d_t> *err;
      arrvec_t<arr_2d_t> *Phi;

      private:

      void cycle_timelevels(int n)
      {
        n = (n+1) % 3;
      }

      void ini_pressure()
      {
	// dt/2 * (Prs-Prs_amb) / rho
	(*Phi)[n-2](this->i^this->halo, this->j^this->halo) = real_t(0);
	(*Phi)[n-1](this->i^this->halo, this->j^this->halo) = real_t(0);
	(*Phi)[n](this->i^this->halo, this->j^this->halo) = real_t(0);
      }

      void pressure_solver_update(real_t dt)
      {
	
        using namespace arakawa_c;
	using formulae::nabla_op::grad;
	using formulae::nabla_op::div;

	real_t beta = .25;  //TODO
        real_t gamma = 1.;  //TODO
	real_t rho = 1.;    //TODO    

        real_t max_err = .0001;

	int halo = this->halo;
	rng_t i = this->i;
	rng_t j = this->j;

	tmp_u = this->psi(u);
	tmp_w = this->psi(w);

    std::cerr<<"--------------------------------------------------------------"<<std::endl;
	//pseudo-time loop
	real_t error = 1.;
	while (error > max_err)
	{
	  this->xchng((*Phi)[n],   i^halo, j^halo);
	  this->xchng((*Phi)[n-1],   i^halo, j^halo);  //some better initial guess needed
	  this->xchng((*Phi)[n-2],   i^halo, j^halo);
	  this->xchng(tmp_u, i^halo, j^halo);
	  this->xchng(tmp_w, i^halo, j^halo);

          //err[n-2]
	  tmp_x(i, j) = rho * tmp_u(i, j) - grad<0>((*Phi)[n-2](i^halo, j^halo), i, j, real_t(1));
	  tmp_z(i, j) = rho * tmp_w(i, j) - grad<1>((*Phi)[n-2](i^halo, j^halo), j, i, real_t(1));
     
	  this->xchng(tmp_x, i^halo, j^halo);
	  this->xchng(tmp_z, i^halo, j^halo);

          (*err)[n-2](i, j) = 1./ rho * div(tmp_x(i^halo,j^halo), tmp_z(i^halo, j^halo), i, j, real_t(1), real_t(1)); //error

          //err[n-1]
	  tmp_x(i, j) = rho * tmp_u(i, j) - grad<0>((*Phi)[n-1](i^halo, j^halo), i, j, real_t(1));
	  tmp_z(i, j) = rho * tmp_w(i, j) - grad<1>((*Phi)[n-1](i^halo, j^halo), j, i, real_t(1));
     
	  this->xchng(tmp_x, i^halo, j^halo);
	  this->xchng(tmp_z, i^halo, j^halo);

          (*err)[n-1](i, j) = 1./ rho * div(tmp_x(i^halo,j^halo), tmp_z(i^halo, j^halo), i, j, real_t(1), real_t(1)); //error

          this->xchng((*err)[n-1], i^halo, j^halo);

          tmp_e1(i, j) = grad<0>((*err)[n-1](i^halo, j^halo), i, j, real_t(1));
          tmp_e2(i, j) = grad<1>((*err)[n-1](i^halo, j^halo), j, i, real_t(1));
          this->xchng(tmp_e1, i^halo, j^halo);
          this->xchng(tmp_e2, i^halo, j^halo);

          lap_err(i,j) = div(tmp_e1(i^halo,j^halo), tmp_e2(i^halo, j^halo), i, j, real_t(1), real_t(1)); //laplasjan(error)

//          tmp_e1(i,j) = err(i,j)*lap_err(i,j);
//          tmp_e2(i,j) = lap_err(i,j)*lap_err(i,j);
//          beta = - blitz::sum(tmp_e1(i,j))/blitz::sum(tmp_e2(i,j));

          (*Phi)[n](i, j) = gamma * ((*Phi)[n-1](i,j) - (*Phi)[n-2](i,j)) + (*Phi)[n-2](i,j) + beta * (*err)[n-1](i, j);
          (*err)[n](i, j) = gamma * ((*err)[n-1](i,j) - (*err)[n-2](i,j)) + (*err)[n-1](i,j) + beta * lap_err(i,j); 
                                                                                                      //lap(err)[n-1]
          error = std::max(std::abs(max((*err)[n])), std::abs(min((*err)[n])));
std::cerr<<"wip !!! "<<std::endl;
          iters++;
          if (error > max_err) cycle_timelevels(n);
	}
	//end of pseudo_time loop
	this->xchng((*Phi)[n], i^halo, j^halo);

	tmp_u(i, j) -= grad<0>((*Phi)[n](i^halo, j^halo), i, j, real_t(1));
	tmp_w(i, j) -= grad<1>((*Phi)[n](i^halo, j^halo), j, i, real_t(1));

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

      struct params_t : parent_t::params_t { };

      // ctor
      pressure_crk2(
	mem_t &mem,
        typename parent_t::bc_p &bcxl,
        typename parent_t::bc_p &bcxr,
        typename parent_t::bc_p &bcyl,
        typename parent_t::bc_p &bcyr,
	const rng_t &i,
	const rng_t &j,
	const params_t &p
      ) :
	parent_t(mem, bcxl, bcxr, bcyl, bcyr, i, j, p),
        n(2),
        // (i, j)
        lap_err(mem.tmp[std::string(__FILE__)][0][0]),
        // (i^hlo, j^hlo))
	tmp_x(mem.tmp[std::string(__FILE__)][0][1]),
	tmp_z(mem.tmp[std::string(__FILE__)][0][2]),
	tmp_u(mem.tmp[std::string(__FILE__)][0][3]),
	tmp_w(mem.tmp[std::string(__FILE__)][0][4]),
	tmp_e1(mem.tmp[std::string(__FILE__)][0][5]),
	tmp_e2(mem.tmp[std::string(__FILE__)][0][6])
      {
        // time levels
        err = &mem.tmp[std::string(__FILE__)][1];
        Phi = &mem.tmp[std::string(__FILE__)][2];
      }

      static void alloc(mem_t &mem, const int nx, const int ny)
      {
        parent_t::alloc(mem, nx, ny);
          
        const std::string file(__FILE__);
        const rng_t i(0, nx-1), j(0, ny-1);
        const int halo = parent_t::halo;
        
        // temporary fields
        mem.tmp[file].push_back(new arrvec_t<arr_2d_t>());
        for (int n=0; n < 1; ++n)   // lap_err
           mem.tmp[file].back().push_back(new arr_2d_t(i, j));  
        for (int n=0; n < 6; ++n)   
           mem.tmp[file].back().push_back(new arr_2d_t( i^halo, j^halo ));  
         
        // vector for err[n-2], err[n-1], err[n]
        mem.tmp[file].push_back(new arrvec_t<arr_2d_t>());
        for (int n=0; n < 3; ++n)
          mem.tmp[file].back().push_back(new arr_2d_t(i^halo, j^halo));           

        // vector for Phi[n-2], Phi[n-1], Phi[n]
        mem.tmp[file].push_back(new arrvec_t<arr_2d_t>());
        for (int n=0; n < 3; ++n)
          mem.tmp[file].back().push_back(new arr_2d_t(i^halo, j^halo));           
      }
    }; 
  }; // namespace solvers
}; // namespace advoocat
