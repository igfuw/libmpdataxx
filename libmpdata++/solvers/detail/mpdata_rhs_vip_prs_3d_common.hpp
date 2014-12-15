/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/formulae/nabla_formulae.hpp>
#include <libmpdata++/solvers/mpdata_rhs_vip.hpp> 

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class mpdata_rhs_vip_prs_3d_common : public mpdata_rhs_vip<ct_params_t>
      {
	using parent_t = mpdata_rhs_vip<ct_params_t>;
        using ix = typename ct_params_t::ix;

        public:
        using real_t = typename ct_params_t::real_t;

	protected:

	// member fields
	const typename ct_params_t::real_t prs_tol;
        int iters = 0;

        typename parent_t::arr_t Phi, tmp_u, tmp_v, tmp_w, err, lap_err, lap_tmp1, lap_tmp2, lap_tmp3;

        auto lap(
          typename parent_t::arr_t &arr, 
          const rng_t &i, 
          const rng_t &j, 
          const rng_t &k, 
          const real_t &dx, 
          const real_t &dy,
          const real_t &dz
        ) return_macro(
          this->xchng_sclr(arr, i^this->halo, j^this->halo, k^this->halo);
          lap_tmp1(i, j, k) = formulae::nabla::grad<0>(arr, i, j, k, dx);
          lap_tmp2(i, j, k) = formulae::nabla::grad<1>(arr, j, k, i, dy);
          lap_tmp3(i, j, k) = formulae::nabla::grad<2>(arr, k, i, j, dz);
          this->xchng_sclr(lap_tmp1, i^this->halo, j^this->halo, k^this->halo, /* deriv = */ true);
          this->xchng_sclr(lap_tmp2, i^this->halo, j^this->halo, k^this->halo, /* deriv = */ true);
          this->xchng_sclr(lap_tmp3, i^this->halo, j^this->halo, k^this->halo, /* deriv = */ true);
          ,
          formulae::nabla::div(lap_tmp1, lap_tmp2, lap_tmp3, i, j, k, dx, dy, dz)
        );

	void ini_pressure()
	{ 
	  const int halo = parent_t::halo;
	  // Phi = dt/2 * (Prs-Prs_amb) / rho 
	  Phi(this->i, this->j, this->k) = real_t(0); // ... but assuming zero perturbation at t=0
	  this->xchng_sclr(Phi, this->i^halo, this->j^halo, this->k^halo);
	}

	virtual void pressure_solver_update() = 0;

	void pressure_solver_apply()
	{
	  const rng_t &i = this->i, &j = this->j, &k = this->k;

	  this->state(ix::u)(i, j, k) += tmp_u(i, j, k);
	  this->state(ix::v)(i, j, k) += tmp_v(i, j, k);
	  this->state(ix::w)(i, j, k) += tmp_w(i, j, k);
	}

        void hook_ante_loop(const int nt)
        {
          parent_t::hook_ante_loop(nt);
	  ini_pressure();
 
          // allow pressure_solver_apply at the first time step
          tmp_u(this->i, this->j, this->k) = 0;
          tmp_v(this->i, this->j, this->k) = 0;
          tmp_w(this->i, this->j, this->k) = 0;
        }

        void hook_ante_step()
        {
          parent_t::hook_ante_step(); // velocity extrapolation + forcings
	  pressure_solver_apply(); 
        }
    
        void hook_post_step()
        {
          parent_t::hook_post_step(); // forcings
	  pressure_solver_update();   // intentionally after forcings (pressure solver must be used after all known forcings are applied)
	  pressure_solver_apply();
        }

	public:

	struct rt_params_t : parent_t::rt_params_t 
        { 
          real_t prs_tol;
        };

	// ctor
	mpdata_rhs_vip_prs_3d_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) : 
	  parent_t(args, p),
          prs_tol(p.prs_tol),
           lap_err(args.mem->tmp[__FILE__][0][0]),
             tmp_u(args.mem->tmp[__FILE__][0][1]),
             tmp_v(args.mem->tmp[__FILE__][0][2]),
             tmp_w(args.mem->tmp[__FILE__][0][3]),
               Phi(args.mem->tmp[__FILE__][0][4]),
               err(args.mem->tmp[__FILE__][0][5]),
	  lap_tmp1(args.mem->tmp[__FILE__][0][6]),
	  lap_tmp2(args.mem->tmp[__FILE__][0][7]),
	  lap_tmp3(args.mem->tmp[__FILE__][0][8])
	{} 

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
	{
	  parent_t::alloc(mem, p);
          parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, 9); // (i^hlo,j^hlo)-sized temporary fields
        }
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
