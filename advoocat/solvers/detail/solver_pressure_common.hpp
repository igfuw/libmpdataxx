/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "../../formulae/nabla_formulae.hpp"
#include "solver_velocity_common.hpp"
#include "../solver_inhomo.hpp"

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      template <class solver_t, int u, int w>
      class pressure_solver_common : public solver_velocity_common<inhomo_solver<solver_t, strang>, u, w>
      {
	protected:

	using parent_t = solver_velocity_common<inhomo_solver<solver_t, solvers::strang>, u, w>;
	typedef typename parent_t::real_t real_t;

	// member fields
	const real_t tol;
        int iters = 0;

        typename parent_t::arr_t Phi, tmp_u, tmp_w, err, lap_err, lap_tmp1, lap_tmp2;

        auto lap(typename parent_t::arr_t &arr, rng_t i, rng_t j, real_t dx, real_t dy) 
        return_macro(
          this->xchng(arr, i^this->halo, j^this->halo);
          lap_tmp1(i, j) = formulae::nabla::grad<0>(arr, i, j, dx);
          lap_tmp2(i, j) = formulae::nabla::grad<1>(arr, j, i, dy);
          this->xchng(lap_tmp1, i^this->halo, j^this->halo);
          this->xchng(lap_tmp2, i^this->halo, j^this->halo);
          ,
          formulae::nabla::div(lap_tmp1, lap_tmp2, i, j, dx, dy)
        );

	void ini_pressure()
	{ 
	  const int halo = parent_t::halo;
	  // Phi = dt/2 * (Prs-Prs_amb) / rho 
	  Phi(this->i, this->j) = real_t(0); // ... but assuming zero perturbation at t=0
	  this->xchng(Phi, this->i^halo, this->j^halo);
	}

	virtual void pressure_solver_update() = 0;

	void pressure_solver_apply()
	{
	  const rng_t &i = this->i, &j = this->j;

	  this->state(u)(i,j) += tmp_u(i,j);
	  this->state(w)(i,j) += tmp_w(i,j);
	}

	public:

        void hook_ante_loop()
        {
          parent_t::hook_ante_loop();
	  ini_pressure();
 
          // allow pressure_solver_apply at the first time step
          tmp_u(this->i, this->j) = 0;
          tmp_w(this->i, this->j) = 0;
        }

        void hook_ante_step()
        {
          parent_t::hook_ante_step(); // velocity extrapolation + forcings
	  pressure_solver_apply(); 
        }
    
        void hook_post_step()
        {
          parent_t::hook_post_step(); // forcings
	  pressure_solver_update(); // intentionally after forcings
	  pressure_solver_apply();
        }

        void hook_post_loop()
        {
          parent_t::hook_post_loop();
std::cerr<<"number of pseudo time iterations "<<iters<<std::endl;
	}

	struct params_t : parent_t::params_t 
        { 
          real_t tol;
        };

	// ctor
	pressure_solver_common(
	  typename parent_t::ctor_args_t args,
	  const params_t &p
	) : 
	  parent_t(args, p),
          tol(p.tol),
          lap_err(args.mem->tmp[__FILE__][0][0]),
          tmp_u(args.mem->tmp[__FILE__][0][1]),
          tmp_w(args.mem->tmp[__FILE__][0][2]),
          Phi(args.mem->tmp[__FILE__][0][3]),
          err(args.mem->tmp[__FILE__][0][4]),
	  lap_tmp1(args.mem->tmp[__FILE__][0][5]),
	  lap_tmp2(args.mem->tmp[__FILE__][0][6])
	{} 

	static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
	{
	  parent_t::alloc(mem, nx, ny);

	  const rng_t i(0, nx-1), j(0, ny-1);
	  const int halo = parent_t::halo;

          // (i^hlo,j^hlo)-sized temporary fields
	  mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
	  for (int n=0; n < 7; ++n)
	    mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t( i^halo, j^halo ));
        }
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
