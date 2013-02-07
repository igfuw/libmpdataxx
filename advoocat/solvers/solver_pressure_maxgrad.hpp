/** 
  * @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once
#include "solver_pressure_common.hpp"
#include "../formulae/nabla_formulae.hpp" //gradient, diveregnce

namespace advoocat
{
  namespace solvers
  {
    template <class inhomo_solver_t, int u, int w, int tht>
    class pressure_maxgrad : public detail::pressure_solver_common<inhomo_solver_t, u, w, tht>
    {
      public:

      using parent_t = detail::pressure_solver_common<inhomo_solver_t, u, w, tht>;
      typedef typename parent_t::mem_t mem_t;
      typedef typename parent_t::real_t real_t;
      using arr_2d_t = typename mem_t::arr_t;

      int iters = 0;

      arr_2d_t Phi; 
      //TODO probably don't need those
      arr_2d_t tmp_u, tmp_w, tmp_x, tmp_z;
      arr_2d_t err, lap_err, tmp_e1, tmp_e2;

      private:

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
	real_t error = 1.;
	while (error > .0001)
	{
	  this->xchng(Phi,   i^halo, j^halo);
	  this->xchng(tmp_u, i^halo, j^halo);
	  this->xchng(tmp_w, i^halo, j^halo);

	  tmp_x(i, j) = rho * tmp_u(i, j) - grad<0>(Phi(i^halo, j^halo), i, j, real_t(1));
	  tmp_z(i, j) = rho * tmp_w(i, j) - grad<1>(Phi(i^halo, j^halo), j, i, real_t(1));
     
	  this->xchng(tmp_x, i^halo, j^halo);
	  this->xchng(tmp_z, i^halo, j^halo);

          err(i, j) = 1./ rho * div(tmp_x(i^halo,j^halo), tmp_z(i^halo, j^halo), i, j, real_t(1), real_t(1)); //error

          this->xchng(err, i^halo, j^halo);

          tmp_e1(i, j) = grad<0>(err(i^halo, j^halo), i, j, real_t(1));
          tmp_e2(i, j) = grad<1>(err(i^halo, j^halo), j, i, real_t(1));
          this->xchng(tmp_e1, i^halo, j^halo);
          this->xchng(tmp_e2, i^halo, j^halo);

          lap_err(i,j) = div(tmp_e1(i^halo,j^halo), tmp_e2(i^halo, j^halo), i, j, real_t(1), real_t(1)); //laplasjan(error)

          tmp_e1(i,j) = err(i,j)*lap_err(i,j);
          tmp_e2(i,j) = lap_err(i,j)*lap_err(i,j);
          beta = - blitz::sum(tmp_e1(i,j))/blitz::sum(tmp_e2(i,j));

          Phi(i, j) -= beta * err(i, j);

          error = std::max(std::abs(max(err)), std::abs(min(err)));
          iters++;
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

      struct params_t : parent_t::params_t { };

      // ctor
      pressure_maxgrad(
	mem_t &mem,
	const rng_t &i,
	const rng_t &j,
	const params_t &p
      ) :
	parent_t(mem, i, j, p),
        // (i, j)
        lap_err(mem.tmp[std::string(__FILE__)][0][0]),
        // (i^hlo, j^hlo))
	err(mem.tmp[std::string(__FILE__)][0][1]),
	tmp_x(mem.tmp[std::string(__FILE__)][0][2]),
	tmp_z(mem.tmp[std::string(__FILE__)][0][3]),
	tmp_u(mem.tmp[std::string(__FILE__)][0][4]),
	tmp_w(mem.tmp[std::string(__FILE__)][0][5]),
	Phi(mem.tmp[std::string(__FILE__)][0][6]),
	tmp_e1(mem.tmp[std::string(__FILE__)][0][7]),
	tmp_e2(mem.tmp[std::string(__FILE__)][0][8])
      {}

      static void alloctmp(
        std::unordered_map<std::string, boost::ptr_vector<arrvec_t<arr_2d_t>>> &tmp,  
        const int nx, 
        const int ny
      )
      {
        parent_t::alloctmp(tmp, nx, ny);

        const rng_t i(0, nx-1), j(0, ny-1);
        const int hlo = 1; // TODO!!!

        // temporary fields
        tmp[std::string(__FILE__)].push_back(new arrvec_t<arr_2d_t>());
        {
          for (int n=0; n < 1; ++n) 
            tmp[std::string(__FILE__)].back().push_back(new arr_2d_t(i, j)); 
          for (int n=0; n < 8; ++n) 
            tmp[std::string(__FILE__)].back().push_back(new arr_2d_t( i^hlo, j^hlo )); 
        }
      }
    }; 
  }; // namespace solvers
}; // namespace advoocat
