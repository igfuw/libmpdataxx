/** 
  * @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  * @brief generalized conjugate residual pressure solver 
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
      template <class ct_params_t, int k_iters, int minhalo>
      class mpdata_rhs_vip_prs_gcrk : public detail::mpdata_rhs_vip_prs_common<ct_params_t, minhalo>
      {
        public:

	using real_t = typename ct_params_t::real_t;

        private:

	using parent_t = detail::mpdata_rhs_vip_prs_common<ct_params_t, minhalo>;
        using ix = typename ct_params_t::ix;

	real_t beta;
        std::vector<real_t> alpha, tmp_den;
	typename parent_t::arr_t lap_err;
	arrvec_t<typename parent_t::arr_t> p_err, lap_p_err;
	
        void pressure_solver_loop_init(bool simple) final
        {
	  p_err[0](this->ijk) = this->err(this->ijk);
	  lap_p_err[0](this->ijk) = this->lap(p_err[0], this->ijk, this->dijk, false, simple);
        }

        void pressure_solver_loop_body(bool simple) final
        {
          for (int v = 0; v < k_iters; ++v)
          {
            tmp_den[v] = this->prs_sum(lap_p_err[v], lap_p_err[v], this->ijk);
            if (tmp_den[v] != 0) beta = - this->prs_sum(this->err, lap_p_err[v], this->ijk) / tmp_den[v];
            this->Phi(this->ijk) += beta * p_err[v](this->ijk);
            this->err(this->ijk) += beta * lap_p_err[v](this->ijk);

            real_t error = std::max(
              std::abs(this->mem->max(this->rank, this->err(this->ijk))), 
              std::abs(this->mem->min(this->rank, this->err(this->ijk)))
            );

            if (error <= this->err_tol) this->converged = true;

            lap_err(this->ijk) = this->lap(this->err, this->ijk, this->dijk, false, simple);

            for (int l = 0; l <= v; ++l)
            {
              if (tmp_den[l] != 0) 
                alpha[l] = - this->prs_sum(lap_err, lap_p_err[l], this->ijk) / tmp_den[l];
            }
            
            if (v < (k_iters - 1))
            {
              p_err[v + 1](this->ijk) = this->err(this->ijk);  
              lap_p_err[v + 1](this->ijk) = lap_err(this->ijk);
              
              for (int l = 0; l <= v; ++l)
              {
                p_err[v + 1](this->ijk) += alpha[l] * p_err[l](this->ijk);
                lap_p_err[v + 1](this->ijk) += alpha[l] * lap_p_err[l](this->ijk);
              }

            }
            else
            {
              p_err[0](this->ijk) = this->err(this->ijk) + alpha[0] * p_err[0](this->ijk);  
              lap_p_err[0](this->ijk) = lap_err(this->ijk) + alpha[0] * lap_p_err[0](this->ijk);
              for (int l = 1; l <= v; ++l)
              {
                p_err[0](this->ijk) += alpha[l] * p_err[l](this->ijk);
                lap_p_err[0](this->ijk) += alpha[l] * lap_p_err[l](this->ijk);
              }
            }
          }
        }

	public:

	struct rt_params_t : parent_t::rt_params_t { };

	// ctor
	mpdata_rhs_vip_prs_gcrk(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
          beta(.25),
          alpha(k_iters, 1.),
          tmp_den(k_iters, 1.),
	  lap_err(args.mem->tmp[__FILE__][0][0]),
	  lap_p_err(args.mem->tmp[__FILE__][1]),
	      p_err(args.mem->tmp[__FILE__][2])
	{}

	static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
	  parent_t::alloc(mem, n_iters);
	  parent_t::alloc_tmp_sclr(mem, __FILE__, 1);
	  parent_t::alloc_tmp_sclr(mem, __FILE__, k_iters);
	  parent_t::alloc_tmp_sclr(mem, __FILE__, k_iters);
	}
      }; 
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
