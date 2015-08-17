/** 
  * @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  * @brief preconditioned conjugate residual pressure solver 
  *   (for more detailed discussion consult Smolarkiewicz & Szmelter 2011
  *    A Nonhydrostatic Unstructured-Mesh Soundproof Model for Simulation of Internal Gravity Waves
  *    Acta Geophysica)
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
      class mpdata_rhs_vip_prs_2d_pc : public detail::mpdata_rhs_vip_prs_2d_common<ct_params_t>
      {
        public:
	
        using real_t = typename ct_params_t::real_t;

        private:

	using parent_t = detail::mpdata_rhs_vip_prs_2d_common<ct_params_t>;
        using ix = typename ct_params_t::ix;

	const int pc_iters;
	real_t beta, alpha, tmp_den;

	using arr_2d_t = typename parent_t::arr_t;
	arr_2d_t p_err, q_err, lap_p_err, lap_q_err;
	arr_2d_t pcnd_err;

	void precond()  //Richardson scheme
	{
	  const rng_t &i = this->i;
	  const rng_t &j = this->j;

	  //initail q_err for preconditioner
	  q_err(this->ijk) = real_t(0);

	  //initail preconditioner error   
	  this->pcnd_err(this->ijk) = this->lap(this->q_err, i, j, this->di, this->dj) - this->err(this->ijk);
	    //TODO does it change with non_const density?
	  
	  assert(pc_iters >= 0 && pc_iters < 10 && "params.pc_iters not specified?");
	  for (int it=0; it<=pc_iters; it++)
	  {
	    q_err(this->ijk)    += real_t(.25) * pcnd_err(this->ijk);
	    pcnd_err(this->ijk) += real_t(.25) * this->lap(this->pcnd_err, i, j, this->di, this->dj);
	  }
	}

        void pressure_solver_loop_init()
        {
	  precond();
	  p_err(this->ijk) = q_err(this->ijk);
	  this->lap_p_err(this->ijk) = this->lap(this->p_err, this->i, this->j, this->di, this->dj);
        }

        void pressure_solver_loop_body()
        {
          tmp_den = this->prs_sum(lap_p_err, lap_p_err, this->i, this->j);
          if (tmp_den != 0) beta = -this->prs_sum(this->err, lap_p_err, this->i, this->j) / tmp_den;
 
          this->Phi(this->ijk) += beta * p_err(this->ijk);
          this->err(this->ijk) += beta * lap_p_err(this->ijk);

          real_t error = std::max(
            std::abs(this->mem->max(this->rank, this->err(this->ijk))), 
            std::abs(this->mem->min(this->rank, this->err(this->ijk)))
          );

          if (error <= this->err_tol) this->converged = true;

          precond();

          this->lap_q_err(this->ijk) = this->lap(this->q_err, this->i, this->j, this->di, this->dj);

          if (tmp_den != 0) alpha = -this->prs_sum(lap_q_err, lap_p_err, this->i, this->j) / tmp_den;

          p_err(this->ijk) *= alpha;
          p_err(this->ijk) += q_err(this->ijk);  
 
          lap_p_err(this->ijk) *= alpha;
          lap_p_err(this->ijk) += lap_q_err(this->ijk);
        }

	public:

	struct rt_params_t : parent_t::rt_params_t { int pc_iters; };

	// ctor
	mpdata_rhs_vip_prs_2d_pc(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
	  parent_t(args, p),
	  pc_iters(p.pc_iters),
          beta(.25),
          alpha(1.),
          tmp_den(1.),
	  lap_p_err(args.mem->tmp[__FILE__][0][0]),
	  lap_q_err(args.mem->tmp[__FILE__][0][1]),
	      p_err(args.mem->tmp[__FILE__][0][2]),
	      q_err(args.mem->tmp[__FILE__][0][3]),
	   pcnd_err(args.mem->tmp[__FILE__][0][4])
	{}

	static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
	  parent_t::alloc(mem, n_iters);
	  parent_t::alloc_tmp_sclr(mem, __FILE__, 5);
	}
      }; 
    } // namespcae detail
  } // namespace solvers
} // namespace libmpdataxx
