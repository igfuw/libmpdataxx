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
*/

#pragma once
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_3d_common.hpp>
#include <libmpdata++/formulae/nabla_formulae.hpp> // gradient, diveregnce

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class mpdata_rhs_vip_prs_3d_cr : public detail::mpdata_rhs_vip_prs_3d_common<ct_params_t>
      {
        public:

	using real_t = typename ct_params_t::real_t;

        private:

	using parent_t = detail::mpdata_rhs_vip_prs_3d_common<ct_params_t>;
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
	  rng_t &k = this->k;

	  this->tmp_u(this->ijk) = this->state(ix::u)(this->ijk);
	  this->tmp_v(this->ijk) = this->state(ix::v)(this->ijk);
	  this->tmp_w(this->ijk) = this->state(ix::w)(this->ijk);

	  this->xchng_sclr(this->Phi,   i^halo, j^halo, k^halo);
	  this->xchng_sclr(this->tmp_u, i^halo, j^halo, k^halo);
	  this->xchng_sclr(this->tmp_v, i^halo, j^halo, k^halo);
	  this->xchng_sclr(this->tmp_w, i^halo, j^halo, k^halo);

	  //initail error   
	  this->err(this->ijk) =
	    - 1./ rho * div(rho * this->tmp_u,
                            rho * this->tmp_v,
                            rho * this->tmp_w,
                            i, j, k,
                            this->di, this->dj, this->dk)
	    + this->lap(this->Phi, i, j, k, this->di, this->dj, this->dk);
	    /* + 1./rho * grad(Phi) * grad(rho) */ // should be added if rho is not constant

	  p_err(this->ijk) = this->err(this->ijk);
	  lap_p_err(this->ijk) = this->lap(p_err, i, j, k, this->di, this->dj, this->dk);

	  //pseudo-time loop
	  this->iters = 0;
	  real_t error = 1.;
	  while (error > this->prs_tol)
	  {
	    tmp_den = this->mem->sum(lap_p_err, lap_p_err, i, j, k);
	    if (tmp_den != 0) beta = - this->mem->sum(this->err, lap_p_err, i, j, k) / tmp_den;
	    this->Phi(this->ijk) += beta * p_err(this->ijk);
	    this->err(this->ijk) += beta * lap_p_err(this->ijk);

	    this->lap_err(this->ijk) = this->lap(this->err, i, j, k, this->di, this->dj, this->dk);         

	    if (tmp_den != 0) alpha = - this->mem->sum(this->lap_err, lap_p_err, i, j, k) / tmp_den;          

	    p_err(this->ijk) *= alpha;
	    p_err(this->ijk) += this->err(this->ijk);  
   
	    lap_p_err(this->ijk) *= alpha;
	    lap_p_err(this->ijk) += this->lap_err(this->ijk);
   
	    error = std::max(
	      std::abs(this->mem->max(this->err(this->ijk))), 
	      std::abs(this->mem->min(this->err(this->ijk)))
	    );
	    this->iters++;
	  }
	  //end of pseudo_time loop

  // TODO: record it
  //std::cerr<<"      number of iterations untill convergence: "<<this->iters<<std::endl;
  //std::cerr<<"      error: "<<error<<std::endl;

	  this->xchng_sclr(this->Phi, i^halo, j^halo, k^halo);

	  this->tmp_u(this->ijk) = - grad<0>(this->Phi, i, j, k, this->di);
	  this->tmp_v(this->ijk) = - grad<1>(this->Phi, j, k, i, this->dj);
	  this->tmp_w(this->ijk) = - grad<2>(this->Phi, k, i, j, this->dk);
	}

	public:

	struct rt_params_t : parent_t::rt_params_t { };

	// ctor
	mpdata_rhs_vip_prs_3d_cr(
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
