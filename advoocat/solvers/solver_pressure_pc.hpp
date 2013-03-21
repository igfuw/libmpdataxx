/** 
  * @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  * @brief conjugate residual pressure solver 
  *   (for more detailed discussion consult Smolarkiewicz & Szmelter 2011
  *    A Nonhydrostatic Unstructured-Mesh Soundproof Model for Simulation of Internal Gravity Waves
  *    Acta Geophysica)
  *
  * @section DERIVATION
  * 
  * for introduction see the derivation of minimal residual pressure solver (solver_pressure_mr.hpp)
  * and conjugate residual pressure solver (solver_pressure_cr.hpp)
  * 
  * in preconditioned solver initail equation \f$ \mathcal{L}(\Phi) - R = 0 \f$ 
  * is replaced by \f$ \mathcal{L}^{\prime}(\Phi) - R^{\prime} = 0 \f$
  * so that convergence of variational schmes is accelerated
  *
  * there is no set of rules explaining how to design a good preconditioner
  * 
  * assume that operator \f$ \mathcal{P} \f$ is our preconditioner
  *
  * In principle \f$ \mathcal{P} \f$ can be any linear operator 
  * such that \f$ \mathcal{L} \mathcal{P} ^{-1} \f$ is negative definite.
  * The goal is to augument the initial recurrence so that it converges faster.
  * On the other hand, for the preconditioner to be useful, the convergence of the 
  * auxiliary problem must be sufficiently rapid to overcome the additional effort
  * of inverting \f$ \mathcal{P} \f$
  * 
  * For the preconditioned conjugate residual scheme we start from 
  *
  * \f$ \mathcal{L}(\Phi) - R = 
  *   \frac{\partial^2 \mathcal{P}(\Phi)}{\partial \tau^2} + \frac{1}{T}\frac{\partial \mathcal{P}(\Phi)}{\partial \tau} \f$
  *
  * acting with \f$ \mathcal{P}^{-1} \f$ on the both sides of the above equation leads to the recurrence formulas
  *
  * \f$ \Phi^{n+1} = \Phi^{n} + \beta^{n} \mathcal{P}^{-1}(p^{n})\f$
  *
  * \f$ r^{n+1} = r^{n} + \beta^{n} \mathcal{L} \mathcal{P}^{-1}(p^{n}) \f$
  *
  * \f$ \mathcal{P}^{-1}(p^{n+1}) = \alpha^{n+1} \mathcal{P}^{-1}(p^{n}) + \mathcal{P}^{-1}(r^{n+1})  \f$
  * 
  * redefining \f$ p_{new} = \mathcal{P}^{-1}(p_{old}) \f$ leads to 
  *
  * \f$ p^{n+1} = \alpha^{n+1} p^{n} + q^{n+1} \f$
  *
  * where \f$ q^{n+1} = \mathcal{P}^{-1}(r^{n+1}) \f$
  *
  * and the coefficients
  *
  * \f$ \beta^{n} = - \frac{<r^{n} \mathcal{L}(p^{n})>}{<\mathcal{L}(p^{n}) \mathcal{L}(p^{n})>} \f$
  *
  * \f$ \alpha^{n+1} = - \frac{<\mathcal{L}(q^{n+1}) \mathcal{L}(p^{n})>}{< \mathcal{L}(p^{n}) \mathcal{L}(p^{n})>} \f$
  *
  * The recurrence of \f$ p \f$ implies a recurence for \f$ \mathcal{L}(p) \f$
  *
  * \f$ \mathcal{L}(p^{n+1}) = \mathcal{L}(q^{n+1}) + \alpha^{n+1} \mathcal{L}(p^{n}) \f$
  * 
  * iterations in pseudo-time stop when residual error is smaller than a given value (for example .0001)
  *
  * Added difficulty of preconditioned scheme stems from the need to solve the auxiliary elliptic problem at initialization
  * \f$ p^{0} = \mathcal{P}^{-1}(r^{0}) \f$ and at each iteration \f$ q^{n+1} = \mathcal{P}^{-1}(r^{n+1}) \f$ .
  * In this approach the inversion of \f$ \mathcal{P} \f$ is done using Richardson scheme (minimum residual scheme with \f$ \beta = .25 \f$).
  * 
  * \f$ e = \mathcal{P}^{-1}(r) \;\;\;\;\;\; / \mathcal{P} \f$
  * 
  * \f$ \mathcal{P}(e) - r = 0 \f$
  *
  * \f$ \mathcal{P}(e) - r = \frac{\partial e}{\partial \xi} \f$

  */

#pragma once
#include <advoocat/solvers/detail/solver_pressure_common.hpp>
#include <advoocat/formulae/nabla_formulae.hpp> //gradient, diveregnce

namespace advoocat
{
  namespace solvers
  {
    template <class solver_t, int u, int w>
    class pressure_pc : public detail::pressure_solver_common<solver_t, u, w>
    {
      public:

      using parent_t = detail::pressure_solver_common<solver_t, u, w>;
      typedef typename parent_t::mem_t mem_t;
      typedef typename parent_t::real_t real_t;
      using arr_2d_t = typename parent_t::arr_t;

      arr_2d_t p_err, q_err, lap_p_err, lap_q_err;
      arr_2d_t pcnd_err;   //TODO is it needed?

      private:

      const int pc_iters;

      void precond()  //Richardson scheme
      {
 	using namespace arakawa_c;
	using formulae::nabla::grad;
	
	int halo = this->halo;
	rng_t &i = this->i;
	rng_t &j = this->j;

        //initail q_err for preconditioner
        q_err(i, j) = real_t(0);
        this->xchng(q_err, i^this->halo, j^this->halo);

       //initail preconditioner error   
        this->pcnd_err(i, j) = this->lap(this->q_err, i, j, this->dx, this->dz) - this->err(i, j);
          //TODO does it change with non_const density?
        this->xchng(pcnd_err, i^halo, j^halo);
        
        assert(pc_iters >= 0 && pc_iters < 10 && "params.pc_iters not specified?");
        for (int it=0; it<=pc_iters; it++)
        {
          q_err(i,j)     += real_t(.25) * pcnd_err(i, j);
          pcnd_err(i, j) += real_t(.25) * this->lap(this->pcnd_err, i, j, this->dx, this->dz);

          this->xchng(q_err, i^halo, j^halo);
        }
      }

      void pressure_solver_update()
      {
	using namespace arakawa_c;
	using formulae::nabla::grad;
	using formulae::nabla::div;

	real_t beta = .25;   //TODO
        real_t alpha = 1.;   //TODO
	real_t rho = 1.;     //TODO    
        real_t tmp_den = 1.; //TODO

	int halo = this->halo;
	rng_t &i = this->i;
	rng_t &j = this->j;

	this->tmp_u(i, j) = this->state(u)(i, j);
	this->tmp_w(i, j) = this->state(w)(i, j);

        this->xchng(this->Phi,   i^halo, j^halo);
        this->xchng(this->tmp_u, i^halo, j^halo);
	this->xchng(this->tmp_w, i^halo, j^halo);

        //initail error   
        this->err(i, j) =
          - 1./ rho * div(rho * this->tmp_u, rho * this->tmp_w , i, j, this->dx, this->dz)
          + this->lap(this->Phi, i, j, this->dx, this->dz);
          /* + 1./rho * grad(Phi) * grad(rho) */ // should be added if rho is not constant

        precond();

        p_err(i, j) = q_err(i, j);
        this->xchng(p_err, i^this->halo, j^this->halo);

        this->lap_p_err(i, j) = this->lap(this->p_err, i, j, this->dx, this->dz);
//std::cerr<<"--------------------------------------------------------------"<<std::endl;
	//pseudo-time loop
	real_t error = 1.;
	while (true)
	{
          tmp_den = this->mem->sum(lap_p_err, lap_p_err, i, j);
          if (tmp_den != 0) beta = - this->mem->sum(this->err, lap_p_err, i, j) / tmp_den;
          //else TODO!
 
          this->Phi(i, j) += beta * p_err(i, j);
          this->err(i, j) += beta * lap_p_err(i, j);

          error = std::max(
            std::abs(this->mem->max(this->err(i, j))), 
            std::abs(this->mem->min(this->err(i, j)))
          );

          if (error <= this->tol) break;
//std::cerr<<"error "<<error<<std::endl;

          //TODO exit pseudotime loop here if <err> < error

          precond();

          this->lap_q_err(i, j) = this->lap(this->q_err, i, j, this->dx, this->dz);

          if (tmp_den != 0) alpha = - this->mem->sum(lap_q_err, lap_p_err, i, j) / tmp_den;

          p_err(i, j) *= alpha;
          p_err(i, j) += q_err(i, j);  
 
          lap_p_err(i, j) *= alpha;
          lap_p_err(i, j) += lap_q_err(i, j);
          this->iters++;
        }
	//end of pseudo_time loop

	this->xchng(this->Phi, i^halo, j^halo);

	this->tmp_u(i, j) -= grad<0>(this->Phi, i, j, this->dx);
	this->tmp_w(i, j) -= grad<1>(this->Phi, j, i, this->dz);

	this->tmp_u(i, j) -= this->state(u)(i, j);
	this->tmp_w(i, j) -= this->state(w)(i, j);
      }

      public:

      struct params_t : parent_t::params_t { int pc_iters; };

      // ctor
      pressure_pc(
	typename parent_t::ctor_args_t args,
	const params_t &p
      ) :
	parent_t(args, p),
        pc_iters(p.pc_iters),
        lap_p_err(args.mem->tmp[__FILE__][0][0]), // TODO: parent has unused lap_err
        lap_q_err(args.mem->tmp[__FILE__][0][1]),
	p_err(args.mem->tmp[__FILE__][0][2]),
	q_err(args.mem->tmp[__FILE__][0][3]),
	pcnd_err(args.mem->tmp[__FILE__][0][4])
      {}

      static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
      {
        parent_t::alloc(mem, nx, ny);
        const rng_t i(0, nx-1), j(0, ny-1);

        // temporary fields
        mem->tmp[__FILE__].push_back(new arrvec_t<arr_2d_t>());
	for (int n=0; n < 5; ++n) 
	  mem->tmp[__FILE__].back().push_back(new arr_2d_t(i^parent_t::halo, j^parent_t::halo)); 
      }
    }; 
  }; // namespace solvers
}; // namespace advoocat
