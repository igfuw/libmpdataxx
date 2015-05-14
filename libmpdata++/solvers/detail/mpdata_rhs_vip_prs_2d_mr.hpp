/** 
  * @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  * @brief minimum residual pressure solver 
  *   (for more detailed discussion consult Smolarkiewicz & Margolin 1994 
  *  Appl. Math and Comp. Sci. 
  *  Variational solver for elliptic problems in atmospheric flows)
  *
  * @section DERIVATION
  * 
  * equations are solved for pressure perturbation (in reference to inertial ambient state)
  * \f$ \Phi = \frac{p-p_e}{\bar{\rho}} \f$
  *
  * where: 
  *
  * \f$ p_e \f$ is pressure of the inertial ambient state (a state that satisfies the equations)
  *
  * \f$ \bar{\rho} \f$ is the reference state density of the fluid
  * 
  * from the continuity equation (applied after the first half-step of advection scheme)
  * 
  * \f$ -\frac{1}{\rho} \nabla \cdot (\rho (\hat{u} - \frac{\triangle t}{2} \nabla \Phi)) = 0 \f$
  *
  * where:
  *
  * \f$ \hat{u} \f$ is the velocity after the first half-step of advection (after applying all the known forcing terms)
  * 
  * multiplying by -1 is to get the diffusion equation
  * 
  * to derive the iterative solver we need to augment continuity equation with a pseudo-time (\f$ \tau \f$) dependence
  * 
  * \f$ -\frac{1}{\rho} \nabla \cdot (\rho (\hat{u} - \frac{\triangle t}{2} \nabla \Phi)) 
      = \frac{\partial \Phi}{\partial \tau} \f$
  *
  * discretizing in pseudo-time (with increment \f$ \beta \f$) leads to
  *
  * \f$ \Phi^{n+1} = \Phi^{n} + \beta (-\frac{1}{\rho} \nabla \cdot (\rho(\hat{u} - \nabla \Phi^{n})) ) \f$
  *
  * where term \f$ -\frac{1}{\rho} \nabla \cdot (\rho(\hat{u} - \nabla \Phi^{n})) \f$ can be denoted as \f$ r^{n} \f$
  * 
  * \f$ r^{n} \f$ is the residual error that needs to be minimized by the pseudo-time iterations within the pressure solver 
  *
  * \f$ \Phi^{n+1} = \Phi^{n} + \beta r^{n} \f$
  *
  * from the definition of \f$ r^{n} \f$ we can write that
  *
  * \f$ r^{n+1} = -\frac{1}{\rho} \nabla \cdot (\rho(\hat{u} - \nabla \Phi^{n+1})) \f$
  *
  * \f$ r^{n+1} = -\frac{1}{\rho} \nabla \cdot (\rho(\hat{u} - \nabla(\Phi^{n} + \beta r^{n}))) \f$
  *
  * \f$ r^{n+1} = r^{n} + \beta \Delta (r^{n}) \f$
  *
  * equations for \f$ \Phi^{n+1} \f$ and \f$ r^{n+1} \f$ form the recurrence formula used by the solver. 
  * To assure convergence pseudo-time increment for solver iterations has to be chosen. 
  * \f$ \beta = const = .25 \f$ assures convergence for every case (this is Richardson scheme).
  *
  * To assure the fastest possible convergence in a given case beta (pseudo-time increment) 
  * has to be chosen in a way that minimizes \f$ r^{n+1} \f$
  * 
  * \f$ <(r^{n+1})^2> \to min \;\;\;\;\;\; / \frac{\partial}{\partial \beta} \;\;\;\;\;\; 
  * \Rightarrow\;\;\;\;\;\; <2(r^{n+1}) \Delta(r^{n})> = 0 \f$
  *
  * where \f$ <r> \f$  denotes the norm of \f$ r \f$ (meaning the sum of \f$ r \f$ over all grid points)
  *
  * substituting \f$ r^{n+1} \f$ by it's recurrence formula and further rearranging leads to the formula for pseudo-time step
  *
  * \f$ \beta = - \frac{r^{n} \Delta r^{n}}{\Delta r^{n} \Delta r^{n}} \f$
  *
  *
  * pseudo-time iterations stop when residual error is smaller than a given value (for example .0001)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_2d_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class mpdata_rhs_vip_prs_2d_mr : public mpdata_rhs_vip_prs_2d_common<ct_params_t>
      {
        public:

        using real_t = typename ct_params_t::real_t;

        private:

	using parent_t = mpdata_rhs_vip_prs_2d_common<ct_params_t>;
        using ix = typename ct_params_t::ix;

	real_t beta, tmp_den;
	typename parent_t::arr_t lap_err;

        void pressure_solver_loop_init() {}

        void pressure_solver_loop_body()
        {
          this->lap_err(this->i, this->j) = this->lap(this->err, this->i, this->j, this->di, this->dj); 

          tmp_den = this->mem->sum(this->lap_err, this->lap_err, this->i, this->j);
          if (tmp_den != 0) beta = - this->mem->sum(this->err, this->lap_err, this->i, this->j) / tmp_den;

          this->Phi(this->ijk) += beta * this->err(this->ijk);
          this->err(this->ijk) += beta * this->lap_err(this->ijk);

          real_t error = std::max(
            std::abs(this->mem->max(this->rank, this->err(this->ijk))), 
            std::abs(this->mem->min(this->rank, this->err(this->ijk)))
          );

          if (error <= this->prs_tol) this->converged = true;
        }

        public:

	struct rt_params_t : parent_t::rt_params_t { };

	// ctor
	mpdata_rhs_vip_prs_2d_mr(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
          beta(.25),
          tmp_den(1.),
	  lap_err(args.mem->tmp[__FILE__][0][0])
	{}

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
	{
	  parent_t::alloc(mem, p);
	  parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, 1);
	}
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
