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
      class mpdata_rhs_vip_prs_common : public mpdata_rhs_vip<ct_params_t>
      {
	using parent_t = mpdata_rhs_vip<ct_params_t>;
        using ix = typename ct_params_t::ix;
	using ijk_t = decltype(mpdata_rhs_vip_prs_common<ct_params_t>::ijk);

        public:
        using real_t = typename ct_params_t::real_t;
        using arr_t = typename parent_t::arr_t;

	protected:

	// member fields
	const real_t prs_tol, err_tol;
        int iters = 0;
        bool converged = false;

        arr_t Phi, err;
        arrvec_t<arr_t> &tmp_uvw, &lap_tmp;

        real_t prs_sum(const arr_t &arr, const ijk_t &ijk)
        {
          return this->mem->sum(arr, ijk, ct_params_t::prs_khn);
        }

        real_t prs_sum(const arr_t &arr1, const arr_t &arr2, const ijk_t &ijk)
        {
          return this->mem->sum(arr1, arr2, ijk, ct_params_t::prs_khn);
        }

        auto lap(
          arr_t &arr, 
          const ijk_t &ijk, 
          const std::array<real_t, parent_t::n_dims>& dijk, 
          bool err_init, // if true then subtract initial state for error calculation
          bool simple // if true do not normalize gradients (simple laplacian)
        ) return_macro(
          this->xchng_pres(arr, ijk);
          formulae::nabla::calc_grad<parent_t::n_dims>(lap_tmp, arr, ijk, dijk);
          if (err_init)
          {
            for (int d = 0; d < parent_t::n_dims; ++d)
            {
              lap_tmp[d](this->ijk) -= tmp_uvw[d](this->ijk);
            }
          }
          if (this->mem->G)
          {
            for (int d = 0; d < parent_t::n_dims; ++d)
            {
              lap_tmp[d](this->ijk) *= (*this->mem->G)(this->ijk);
            }
          }
          if (!simple) this->normalize_vip(lap_tmp);
          this->set_edges(lap_tmp, this->ijk, err_init ? -1 : 0);
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            this->xchng_pres(lap_tmp[d], ijk);
          }
          ,
          formulae::nabla::div<parent_t::n_dims>(lap_tmp, ijk, dijk)
          / formulae::G<ct_params_t::opts>(*this->mem->G, this->ijk)
        )

	void ini_pressure()
	{ 
	  Phi(this->ijk) = 0;
          int npoints = 1;
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
	    Phi(this->ijk) -= 0.5 * pow2(this->vips()[d](this->ijk));
            npoints *= (this->mem->grid_size[d].last() + 1);
          }
          
          auto Phi_mean = prs_sum(Phi, this->ijk) / npoints;
	  Phi(this->ijk) -= Phi_mean;
	}

	virtual void pressure_solver_loop_init(bool) = 0;
	virtual void pressure_solver_loop_body(bool) = 0;

	void pressure_solver_update(bool simple = false)
        {
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            tmp_uvw[d](this->ijk) = this->vips()[d](this->ijk);
          }

	  //initial error   
          err(this->ijk) = lap(Phi, this->ijk, this->dijk, true, simple);

	  iters = 0;
          converged = false;

          pressure_solver_loop_init(simple);
	  //pseudo-time loop
	  while (!converged)
	  {
            pressure_solver_loop_body(simple);
	    iters++;
          }

	  this->xchng_pres(this->Phi, this->ijk);

          formulae::nabla::calc_grad<parent_t::n_dims>(tmp_uvw, Phi, this->ijk, this->dijk);
        }

	void pressure_solver_apply()
	{
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
	    this->vips()[d](this->ijk) -= tmp_uvw[d](this->ijk);
          }
	}

        void hook_ante_loop(const int nt)
        {
          // save initial edge velocities
          this->save_edges(this->vips(), this->ijk);
	  
          // correct initial velocity
	  Phi(this->ijk) = real_t(0);
	  this->xchng_pres(Phi, this->ijk);

	  pressure_solver_update(true);

          this->xchng_pres(this->Phi, this->ijk);
          formulae::nabla::calc_grad<parent_t::n_dims>(tmp_uvw, Phi, this->ijk, this->dijk);
	  pressure_solver_apply();
          this->set_edges(this->vips(), this->ijk, 1);
	  
          parent_t::hook_ante_loop(nt);

          // potential pressure
          ini_pressure();
 
          // allow pressure_solver_apply at the first time step
	  this->xchng_pres(this->Phi, this->ijk);
          formulae::nabla::calc_grad<parent_t::n_dims>(tmp_uvw, Phi, this->ijk, this->dijk);
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            this->vip_rhs[d](this->ijk) -= tmp_uvw[d](this->ijk);
          }
        }

        void vip_rhs_impl_fnlz()
        {
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            this->vip_rhs[d](this->ijk) = -this->vips()[d](this->ijk);
          }

          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl) this->add_relax();
          pressure_solver_update();   // intentionally after forcings (pressure solver must be used after all known forcings are applied)
          pressure_solver_apply();
          this->normalize_vip(this->vips());
          this->set_edges(this->vips(), this->ijk, 1);

          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            this->vip_rhs[d](this->ijk) += this->vips()[d](this->ijk);
            this->vip_rhs[d](this->ijk) /= (0.5 * this->dt);
          }
        }

	public:

	struct rt_params_t : parent_t::rt_params_t 
        { 
          real_t prs_tol;
        };

	// ctor
	mpdata_rhs_vip_prs_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) : 
	  parent_t(args, p),
          prs_tol(p.prs_tol),
          err_tol(p.prs_tol / this->dt), // make stopping criterion correspond to dimensionless divergence
               Phi(args.mem->tmp[__FILE__][0][0]),
               err(args.mem->tmp[__FILE__][0][1]),
           tmp_uvw(args.mem->tmp[__FILE__][1]),
	   lap_tmp(args.mem->tmp[__FILE__][2])
	{} 

	static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
	  parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // Phi, err
          parent_t::alloc_tmp_sclr(mem, __FILE__, parent_t::n_dims); // tmp_uvw
          parent_t::alloc_tmp_sclr(mem, __FILE__, parent_t::n_dims); // lap_tmp
        }
      }; 
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
