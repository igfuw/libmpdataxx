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
      class mpdata_rhs_vip_prs_2d_cr : public detail::mpdata_rhs_vip_prs_2d_common<ct_params_t>
      {
        public:

	using real_t = typename ct_params_t::real_t;

        private:

	using parent_t = detail::mpdata_rhs_vip_prs_2d_common<ct_params_t>;
        using ix = typename ct_params_t::ix;

	real_t beta, alpha, tmp_den;
	typename parent_t::arr_t lap_err, p_err, lap_p_err;
	
        void pressure_solver_loop_init()
        {
	  p_err(this->ijk) = this->err(this->ijk);
	  lap_p_err(this->ijk) = this->lap(p_err, this->i, this->j, this->di, this->dj);
        }

        void pressure_solver_loop_body()
        {
          tmp_den = this->mem->sum(lap_p_err, lap_p_err, this->i, this->j);
          if (tmp_den != 0) beta = - this->mem->sum(this->err, lap_p_err, this->i, this->j) / tmp_den;
          this->Phi(this->ijk) += beta * p_err(this->ijk);
          this->err(this->ijk) += beta * lap_p_err(this->ijk);

          this->lap_err(this->ijk) = this->lap(this->err, this->i, this->j, this->di, this->dj);         

          if (tmp_den != 0) alpha = - this->mem->sum(this->lap_err, lap_p_err, this->i, this->j) / tmp_den;          

          p_err(this->ijk) *= alpha;
          p_err(this->ijk) += this->err(this->ijk);  
 
          lap_p_err(this->ijk) *= alpha;
          lap_p_err(this->ijk) += this->lap_err(this->ijk);

          real_t error = std::max(
            std::abs(this->mem->max(this->rank, this->err(this->ijk))), 
            std::abs(this->mem->min(this->rank, this->err(this->ijk)))
          );

          if (error <= this->prs_tol) this->converged = true;
        }

	public:

	struct rt_params_t : parent_t::rt_params_t { };

	// ctor
	mpdata_rhs_vip_prs_2d_cr(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
          beta(.25),
          alpha(1.),
          tmp_den(1.),
	  lap_err(args.mem->tmp[__FILE__][0][0]),
	  lap_p_err(args.mem->tmp[__FILE__][0][1]),
	      p_err(args.mem->tmp[__FILE__][0][2])
	{}

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
	{
	  parent_t::alloc(mem, p);
	  parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, 3);
	}
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
