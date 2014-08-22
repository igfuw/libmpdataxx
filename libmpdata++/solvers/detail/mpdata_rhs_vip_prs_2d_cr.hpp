/** 
  * @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  * @brief conjugate residual pressure solver 
  *   (for more detailed discussion consult Smolarkiewicz & Margolin 1994 
  *  Appl. Math and Comp. Sci. 
  *  Variational solver for elliptic problems in atmospheric flows)
  *
  * @section DERIVATION
  * 
  * for introduction see the derivation of minimal residual pressure solver (solver_pressure_mr.hpp)
  * 
  * \f$ -\frac{1}{\bar{\rho}} \nabla \cdot (\bar{\rho} (\hat{u} - \frac{\triangle t}{2} \nabla \Phi)) = 0 \f$
  *
  * above equation can be written as \f$ \mathcal{L}(\Phi) - R = 0 \f$ 
  * 
  * where \f$ \mathcal{L}() \f$ in theory may be any linear semidefinite operator
  * (this scheme doesn't require the operator to be self-adjoint)
  *
  * (for a concise discussion of the needed assumtions for operator \f$ \mathcal{L}() \f$ consult Smolarkiewicz & Margolin 1994)
  *
  * in this case:
  *
  * \f$ \mathcal{L}() = \Delta() \f$
  *
  * \f$ R = - \frac{1}{\rho} \nabla \cdot {\rho} \hat{u} \f$ 
  *
  * to obtain faster convergence (than minimum residual scheme) we start from dampened wave equation (instead of diffusion equation)
  * 
  * \f$ \mathcal{L}(\Phi) - R = 0 \;\;\;\;\;\;   \Rightarrow \;\;\;\;\;\;\;\;
  * \mathcal{L}(\Phi) - R = \frac{\partial^2 \Phi}{\partial \tau^2} + \frac{1}{T}\frac{\partial \Phi}{\partial \tau} \f$ 
  *
  * using centered differencing for the second derivative and one-sided differencing or the first derivative 
  * leads to three term recurrence formula
  *
  * \f$ \Phi^{n+1} = \gamma \Phi^{n} + (1-\gamma)\Phi^{n-1} + \beta (\mathcal{L}(\Phi^{n}) - R)  \f$
  *
  * where
  * 
  * \f$ \gamma = \frac{2+\frac{\triangle \tau}{T}}{1+\frac{\triangle \tau}{T}}  \f$
  *
  * \f$ \beta = \frac{\triangle \tau^2}{1+\frac{\triangle \tau}{T}}  \f$
  * 
  * which after rearranging leads to
  *
  * \f$ \Phi^{n+1} = \Phi^{n} + \beta ^{n} (\alpha ^{n} p^{n-1} + r^{n}) \f$
  * 
  * where:
  *
  * \f$ \alpha ^{n} = \frac{(\gamma ^{n} -1) \beta ^{n-1}}{\beta ^{n}} \f$
  *
  * \f$ p^{n} = \frac{\Phi^{n+1}-\Phi ^{n}}{\beta ^{n} } \f$
  *
  * \f$ r^n = \mathcal{L}(\Phi ^{n}) -R \f$
  *
  * this leads to recurrence algorithm for \f$ \Phi \f$, residual error (\f$ r \f$) and directional error (\f$ p \f$)
  * 
  * \f$ \Phi ^{n+1} = \Phi ^{n} + \beta ^{n} p^{n} \f$
  *
  * \f$ r^{n+1} = r^{n} + \beta ^{n} \mathcal{L}(p^{n}) \f$
  *
  * \f$ p^{n+1} = r^{n+1} + \alpha ^{n+1} p^{n} \f$
  *
  * with the coefficients
  *
  * \f$ \beta ^{n}= - \frac{<r^{n} \mathcal{L}(p^n)>}{<\mathcal{L}(p^n) \mathcal{L}(p^n)>} \f$
  *
  * \f$ \alpha ^{n+1}= -  \frac{<\mathcal{L}(r^{n+1}) \mathcal{L}(p^n)>}{<\mathcal{L}(p^n) \mathcal{L}(p^n)> } \f$
  *
  * The recurrence of \f$ p \f$ implies a recuurence for \f$ \mathcal{L}(p) \f$
  *
  * \f$ \mathcal{L}(p^{n+1}) = \mathcal{L}(r^{n+1}) + \alpha^{n+1} \mathcal{L}(p^{n}) \f$
  * 
  * iterations in pseudo-time stop when residual error is smaller than a given value (for example .0001)

*/

#pragma once
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_2d_common.hpp>
#include <libmpdata++/formulae/nabla_formulae.hpp> // gradient, diveregnce

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class mpdata_rhs_vip_prs_2d_cr : public detail::mpdata_rhs_vip_prs_2d_common<ct_params_t>
      {
	using parent_t = detail::mpdata_rhs_vip_prs_2d_common<ct_params_t>;
	using real_t = typename ct_params_t::real_t;
        using ix = typename ct_params_t::ix;

	typename parent_t::arr_t p_err, lap_p_err;

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

	  this->tmp_u(i, j) = this->state(ix::u)(i, j);
	  this->tmp_w(i, j) = this->state(ix::w)(i, j);

	  //initial error   
          this->err(i, j) = this->err_init(this->Phi, this->tmp_u, this->tmp_w, i, j, this->di, this->dj);
	    /* + 1./rho * grad(Phi) * grad(rho) */ // should be added if rho is not constant

	  p_err(i ,j) = this->err(i, j);
	  lap_p_err(i,j) = this->lap(p_err, i, j, this->di, this->dj);

	  //pseudo-time loop
	  this->iters = 0;
	  real_t error = 1.;
	  while (error > this->prs_tol)
	  {
	    tmp_den = this->mem->sum(lap_p_err, lap_p_err, i, j);
	    if (tmp_den != 0) beta = - this->mem->sum(this->err, lap_p_err, i, j) / tmp_den;
	    this->Phi(i, j) += beta * p_err(i, j);
	    this->err(i, j) += beta * lap_p_err(i, j);

	    this->lap_err(i, j) = this->lap(this->err, i, j, this->di, this->dj);         

	    if (tmp_den != 0) alpha = - this->mem->sum(this->lap_err, lap_p_err, i, j) / tmp_den;          

	    p_err(i, j) *= alpha;
	    p_err(i, j) += this->err(i, j);  
   
	    lap_p_err(i,j) *= alpha;
	    lap_p_err(i,j) += this->lap_err(i,j);
   
	    error = std::max(
	      std::abs(this->mem->max(this->err(i,j))), 
	      std::abs(this->mem->min(this->err(i,j)))
	    );
	    this->iters++;
	  }
	  //end of pseudo_time loop

  // TODO: record it
  //std::cerr<<"      number of iterations untill convergence: "<<this->iters<<std::endl;
  //std::cerr<<"      error: "<<error<<std::endl;

	  this->xchng_pres(this->Phi, i^halo, j^halo);

	  this->tmp_u(i, j) = - grad<0>(this->Phi, i, j, this->di);
	  this->tmp_w(i, j) = - grad<1>(this->Phi, j, i, this->dj);
          
          this->set_edges(this->tmp_u, this->tmp_w, this->state(ix::u), this->state(ix::w), i, j);
	}

	public:

	struct rt_params_t : parent_t::rt_params_t { };

	// ctor
	mpdata_rhs_vip_prs_2d_cr(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
	  lap_p_err(args.mem->tmp[__FILE__][0][0]),
	      p_err(args.mem->tmp[__FILE__][0][1])
	{}

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
	{
	  parent_t::alloc(mem, p);
	  parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, 2);
	}
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
