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
*/

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t, int minhalo>
      class mpdata_rhs_vip_prs_mr : public mpdata_rhs_vip_prs_common<ct_params_t, minhalo>
      {
        public:

        using real_t = typename ct_params_t::real_t;

        private:

	using parent_t = mpdata_rhs_vip_prs_common<ct_params_t, minhalo>;
        using ix = typename ct_params_t::ix;

	real_t beta, tmp_den;
	typename parent_t::arr_t lap_err;

        void pressure_solver_loop_init(bool simple) final {}

        void pressure_solver_loop_body(bool simple) final
        {
          this->lap_err(this->ijk) = this->lap(this->err, this->ijk, this->dijk, false, simple);

          tmp_den = this->prs_sum(this->lap_err, this->lap_err, this->ijk);
          if (tmp_den != 0) beta = - this->prs_sum(this->err, this->lap_err, this->ijk) / tmp_den;

          this->Phi(this->ijk) += beta * this->err(this->ijk);
          this->err(this->ijk) += beta * this->lap_err(this->ijk);

          real_t error = std::max(
            std::abs(this->mem->max(this->rank, this->err(this->ijk))), 
            std::abs(this->mem->min(this->rank, this->err(this->ijk)))
          );

          if (error <= this->err_tol) this->converged = true;
        }

        public:

	struct rt_params_t : parent_t::rt_params_t { };

	// ctor
	mpdata_rhs_vip_prs_mr(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
          beta(.25),
          tmp_den(1.),
	  lap_err(args.mem->tmp[__FILE__][0][0])
	{}

	static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
	  parent_t::alloc(mem, n_iters);
	  parent_t::alloc_tmp_sclr(mem, __FILE__, 1);
	}
      }; 
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
