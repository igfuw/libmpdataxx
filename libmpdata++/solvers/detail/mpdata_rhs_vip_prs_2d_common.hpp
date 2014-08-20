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
      class mpdata_rhs_vip_prs_2d_common : public mpdata_rhs_vip<ct_params_t>
      {
	using parent_t = mpdata_rhs_vip<ct_params_t>;
        using ix = typename ct_params_t::ix;

        public:
        using real_t = typename ct_params_t::real_t;

	protected:

	// member fields
	const typename ct_params_t::real_t prs_tol;
        int iters = 0;

        typename parent_t::arr_t Phi, tmp_u, tmp_w, err, lap_err, lap_tmp1, lap_tmp2;

        auto lap(
          typename parent_t::arr_t &arr, 
          const rng_t &i, 
          const rng_t &j, 
          const real_t &dx, 
          const real_t &dy
        ) return_macro(
          this->xchng_pres(arr, i, j);
          lap_tmp1(i, j) = formulae::nabla::grad<0>(arr, i, j, dx);
          lap_tmp2(i, j) = formulae::nabla::grad<1>(arr, j, i, dy);
          this->set_edges(lap_tmp1, lap_tmp2, i, j);
          this->xchng_pres(lap_tmp1, i, j);
          this->xchng_pres(lap_tmp2, i, j);
          ,
          formulae::nabla::div(lap_tmp1, lap_tmp2, i, j, dx, dy)
        );
        
        auto err_init(
          typename parent_t::arr_t &arr, 
          typename parent_t::arr_t &v1, 
          typename parent_t::arr_t &v2, 
          const rng_t &i, 
          const rng_t &j, 
          const real_t &dx, 
          const real_t &dy
        ) return_macro(
          this->xchng_pres(arr, i^this->halo, j^this->halo);
          lap_tmp1(i, j) = formulae::nabla::grad<0>(arr, i, j, dx) - v1(i, j);
          lap_tmp2(i, j) = formulae::nabla::grad<1>(arr, j, i, dy) - v2(i, j);
          this->set_edges(lap_tmp1, lap_tmp2, i, j);
          this->xchng_pres(lap_tmp1, i, j);
          this->xchng_pres(lap_tmp2, i, j);
          ,
          formulae::nabla::div(lap_tmp1, lap_tmp2, i, j, dx, dy)
        );

	void ini_pressure()
	{ 
	  const int halo = parent_t::halo;
	  // Phi = dt/2 * (Prs-Prs_amb) / rho 
	  Phi(this->i, this->j) = real_t(0); // ... but assuming zero perturbation at t=0
	  this->xchng_pres(Phi, this->i^halo, this->j^halo);
	}
	
        void xchng_pres(typename parent_t::arr_t arr, rng_t range_i, rng_t range_j)
        {
          this->mem->barrier();
          this->bcxl->fill_halos_pres(arr, range_j);
          this->bcxr->fill_halos_pres(arr, range_j);
          this->bcyl->fill_halos_pres(arr, range_i);
          this->bcyr->fill_halos_pres(arr, range_i);
          this->mem->barrier();
        }
        
        void set_edges(typename parent_t::arr_t arr1,
                       typename parent_t::arr_t arr2,
                       rng_t range_i,
                       rng_t range_j)
        {
          this->bcxl->set_edge_pres(arr1, range_j);
          this->bcxr->set_edge_pres(arr1, range_j);
          this->bcyl->set_edge_pres(arr2, range_i);
          this->bcyr->set_edge_pres(arr2, range_i);
          this->mem->barrier();
        }
        
        void set_edges(typename parent_t::arr_t arr1,
                       typename parent_t::arr_t arr2,
                       typename parent_t::arr_t v1,
                       typename parent_t::arr_t v2,
                       rng_t range_i,
                       rng_t range_j)
        {
          this->bcxl->set_edge_pres(arr1, v1, range_j);
          this->bcxr->set_edge_pres(arr1, v1, range_j);
          this->bcyl->set_edge_pres(arr2, v2, range_i);
          this->bcyr->set_edge_pres(arr2, v2, range_i);
          this->mem->barrier();
        }

	virtual void pressure_solver_update() = 0;

	void pressure_solver_apply()
	{
	  const rng_t &i = this->i, &j = this->j;

	  this->state(ix::u)(i,j) += tmp_u(i,j);
	  this->state(ix::w)(i,j) += tmp_w(i,j);
	}

        void hook_ante_loop(const int nt)
        {
          parent_t::hook_ante_loop(nt);
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
	  pressure_solver_update();   // intentionally after forcings (pressure solver must be used after all known forcings are applied)
	  pressure_solver_apply();
        }

	public:

	struct rt_params_t : parent_t::rt_params_t 
        { 
          real_t prs_tol;
        };

	// ctor
	mpdata_rhs_vip_prs_2d_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) : 
	  parent_t(args, p),
          prs_tol(p.prs_tol),
           lap_err(args.mem->tmp[__FILE__][0][0]),
             tmp_u(args.mem->tmp[__FILE__][0][1]),
             tmp_w(args.mem->tmp[__FILE__][0][2]),
               Phi(args.mem->tmp[__FILE__][0][3]),
               err(args.mem->tmp[__FILE__][0][4]),
	  lap_tmp1(args.mem->tmp[__FILE__][0][5]),
	  lap_tmp2(args.mem->tmp[__FILE__][0][6])
	{} 

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
	{
	  parent_t::alloc(mem, p);
          parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, 7); // (i^hlo,j^hlo)-sized temporary fields
        }
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
