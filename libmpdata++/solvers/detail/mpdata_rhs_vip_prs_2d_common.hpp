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
        using arr_t = typename parent_t::arr_t;

	protected:

	// member fields
	const real_t err_tol;
        int iters = 0;
        bool converged = false;

        arr_t Phi, tmp_u, tmp_w, err, lap_tmp1, lap_tmp2;

        real_t prs_sum(const arr_t &arr, const rng_t &i, const rng_t &j)
        {
          return this->mem->sum(arr, i, j, ct_params_t::prs_khn);
        }

        real_t prs_sum(const arr_t &arr1, const arr_t &arr2, const rng_t &i, const rng_t &j)
        {
          return this->mem->sum(arr1, arr2, i, j, ct_params_t::prs_khn);
        }
      
        auto lap(
          arr_t &arr, 
          const rng_t &i, 
          const rng_t &j, 
          const real_t dx, 
          const real_t dy
        ) return_macro(
          this->xchng_pres(arr, i, j);
          lap_tmp1(this->ijk) = formulae::nabla::grad<0>(arr, i, j, dx);
          lap_tmp2(this->ijk) = formulae::nabla::grad<1>(arr, j, i, dy);
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl)
          {
            lap_tmp1(this->ijk) /= (1 + 0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk));
            lap_tmp2(this->ijk) /= (1 + 0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk));
          }
          this->set_edges(lap_tmp1, lap_tmp2, i, j, 0);
          this->xchng_pres(lap_tmp1, i, j);
          this->xchng_pres(lap_tmp2, i, j);
          ,
          formulae::nabla::div(lap_tmp1, lap_tmp2, i, j, dx, dy)
        );
        
        auto err_init(
          arr_t &arr, 
          const arr_t &v1, 
          const arr_t &v2, 
          const rng_t &i, 
          const rng_t &j, 
          const real_t dx, 
          const real_t dy
        ) return_macro(
          this->xchng_pres(arr, i^this->halo, j^this->halo);
          lap_tmp1(this->ijk) = formulae::nabla::grad<0>(arr, i, j, dx) - v1(this->ijk);
          lap_tmp2(this->ijk) = formulae::nabla::grad<1>(arr, j, i, dy) - v2(this->ijk);
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl)
          {
            lap_tmp1(this->ijk) /= (1 + 0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk));
            lap_tmp2(this->ijk) /= (1 + 0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk));
          }
          this->set_edges(lap_tmp1, lap_tmp2, i, j, -1);
          this->xchng_pres(lap_tmp1, i, j);
          this->xchng_pres(lap_tmp2, i, j);
          ,
          formulae::nabla::div(lap_tmp1, lap_tmp2, i, j, dx, dy)
        );

	void ini_pressure()
	{ 
	  const int halo = parent_t::halo;
	  // Phi = dt/2 * (Prs-Prs_amb) / rho 
	  Phi(this->ijk) = real_t(0); // ... but assuming zero perturbation at t=0
	  this->xchng_pres(Phi, this->i^halo, this->j^halo);
	}
	
        void xchng_pres(arr_t &arr, const rng_t &range_i, const rng_t &range_j)
        {
          this->mem->barrier();
          this->bcxl->fill_halos_pres(arr, range_j);
          this->bcxr->fill_halos_pres(arr, range_j);
          this->bcyl->fill_halos_pres(arr, range_i);
          this->bcyr->fill_halos_pres(arr, range_i);
          this->mem->barrier();
        }
        
        void set_edges(arr_t &arr1,
                       arr_t &arr2,
                       const rng_t &range_i,
                       const rng_t &range_j,
                       int sign
        )
        {
          this->bcxl->set_edge_pres(arr1, range_j, sign);
          this->bcxr->set_edge_pres(arr1, range_j, sign);
          this->bcyl->set_edge_pres(arr2, range_i, sign);
          this->bcyr->set_edge_pres(arr2, range_i, sign);
          this->mem->barrier();
        }
        
        void save_edges(arr_t &arr1,
                        arr_t &arr2,
                        const rng_t &range_i,
                        const rng_t &range_j
        )
        {
          this->bcxl->save_edge_vel(arr1, range_j);
          this->bcxr->save_edge_vel(arr1, range_j);
          this->bcyl->save_edge_vel(arr2, range_i);
          this->bcyr->save_edge_vel(arr2, range_i);
        }
        
        virtual void pressure_solver_loop_init() = 0;
        virtual void pressure_solver_loop_body() = 0;

	void pressure_solver_update()
        {
          const auto &i = this->i, &j = this->j;

	  tmp_u(this->ijk) = this->state(ix::vip_i)(this->ijk);
	  tmp_w(this->ijk) = this->state(ix::vip_j)(this->ijk);

	  //initial error   
          err(this->ijk) = err_init(Phi, tmp_u, tmp_w, i, j, this->di, this->dj);

	  iters = 0;
          converged = false;

          pressure_solver_loop_init();
	  //pseudo-time loop
	  while (!converged)
	  {
            pressure_solver_loop_body();
	    iters++;
          }

	  xchng_pres(this->Phi, i^this->halo, j^this->halo);

	  using formulae::nabla::grad;
	  tmp_u(this->ijk) = - grad<0>(Phi, i, j, this->di);
	  tmp_w(this->ijk) = - grad<1>(Phi, j, i, this->dj);
        }

	void pressure_solver_apply()
	{
	  this->state(ix::vip_i)(this->ijk) += tmp_u(this->ijk);
	  this->state(ix::vip_j)(this->ijk) += tmp_w(this->ijk);
          set_edges(this->state(ix::vip_i), this->state(ix::vip_j), this->i, this->j, 1);
	}

        void hook_ante_loop(const real_t tshift)
        {
          parent_t::hook_ante_loop(tshift);
	  ini_pressure();
          
          // save initial edge velocities
          save_edges(this->state(ix::vip_i), this->state(ix::vip_j), this->i, this->j);
 
          // allow pressure_solver_apply at the first time step
          tmp_u(this->ijk) = 0;
          tmp_w(this->ijk) = 0;
        }

        void hook_ante_step()
        {
          parent_t::hook_ante_step(); // velocity extrapolation + forcings
	  pressure_solver_apply(); 
        }
    
        void hook_post_step()
        {
          // intentionally calling "great-grandparent" i.e. mpdata_rhs hook
          parent_t::parent_t::parent_t::hook_post_step(); // forcings
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl) parent_t::add_relax();
	  pressure_solver_update();   // intentionally after forcings (pressure solver must be used after all known forcings are applied)
	  pressure_solver_apply();
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl) parent_t::normalize_absorber();
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
          err_tol(p.prs_tol / this->dt), // make stopping criterion correspond to dimensionless divergence
             tmp_u(args.mem->tmp[__FILE__][0][0]),
             tmp_w(args.mem->tmp[__FILE__][0][1]),
               Phi(args.mem->tmp[__FILE__][0][2]),
               err(args.mem->tmp[__FILE__][0][3]),
	  lap_tmp1(args.mem->tmp[__FILE__][0][4]),
	  lap_tmp2(args.mem->tmp[__FILE__][0][5])
	{} 

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
	{
	  parent_t::alloc(mem, p);
          parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, 6); // (i^hlo,j^hlo)-sized temporary fields
        }
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
