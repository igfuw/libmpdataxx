/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    enum sgs_scheme_t
    {
      dns,
      smg,
      tke
    };

    enum stress_diff_t
    {
      normal,
      pade
    };

    // to be specialised
    template <class ct_params_t, class enableif = void>
    class mpdata_rhs_vip_prs_sgs
    {};

    // 2D version
    template <class ct_params_t>
    class mpdata_rhs_vip_prs_sgs<
      ct_params_t,
      typename std::enable_if<ct_params_t::n_dims == 2>::type
    > : public mpdata_rhs_vip_prs<ct_params_t>
    {
      using parent_t = mpdata_rhs_vip_prs<ct_params_t>;

      protected:

      // member fields
      arrvec_t<typename parent_t::arr_t> &tau;
      arrvec_t<typename parent_t::arr_t> &drv;
      arrvec_t<typename parent_t::arr_t> &wrk;
      typename ct_params_t::real_t eta = 0;

      template <int d, class arg_t>
      inline auto pade_helper(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j
      ) return_macro(,
	  x(idxperm::pi<d>(i + 1, j)) + 2 * x(idxperm::pi<d>(i, j)) + x(idxperm::pi<d>(i - 1, j))
      )

      // apply pade scheme to the d-th element of the drv array
      void pade_deriv(int d)
      {
        wrk[0](this->ijk) = drv[d](this->ijk);

        for (int m = 0; m < 3; ++m)
        {
          this->xchng_sclr(wrk[0], this->ijk);
          switch(d)
          {
            case 0 : case 1 :
            {
              wrk[1](this->ijk) = pade_helper<0>(wrk[0], this->i, this->j);
              break;
            }
            case 2 : case 3 :
            {
              wrk[1](this->ijk) = pade_helper<1>(wrk[0], this->j, this->i);
              break;
            }
            default : assert(false);
          }
          wrk[0](this->ijk) += (drv[d](this->ijk) - 0.25 * wrk[1](this->ijk));
        }

        drv[d](this->ijk) = wrk[0](this->ijk);
      }

      void hook_ante_step()
      {
        using namespace formulae::nabla;
        using ix = typename ct_params_t::ix;

        // to avoid writing this->
        const auto &i = this->i, &j = this->j;
        const auto &di = this->di, &dj = this->dj;
        const auto &ijk = this->ijk;

        // TODO: get rid of superfluous barriers
        this->xchng_sclr(this->state(ix::vip_i), ijk);
        this->xchng_sclr(this->state(ix::vip_j), ijk);

        // velocity gradient
        drv[0](ijk) = grad<0>(this->state(ix::vip_i), i, j, di);
        drv[1](ijk) = grad<0>(this->state(ix::vip_j), i, j, di);

        drv[2](ijk) = grad<1>(this->state(ix::vip_i), j, i, dj);
        drv[3](ijk) = grad<1>(this->state(ix::vip_j), j, i, dj);

        if ((stress_diff_t)ct_params_t::stress_diff == pade)
        {
          for (int d = 0; d < 4; ++d)
            pade_deriv(d);
        }

        tau[0](ijk) = drv[0](ijk);
        tau[1](ijk) = 0.5 * (drv[1](ijk) + drv[2](ijk));
        tau[2](ijk) = drv[3](ijk);

        // TODO: get rid of superfluous barriers
        for (auto& t : tau)
          this->xchng_sclr(t, ijk);

        // elements of strain rate tensor divergence
        drv[0](ijk) = grad<0>(tau[0], i, j, di);
        drv[1](ijk) = grad<0>(tau[1], i, j, di);

        drv[2](ijk) = grad<1>(tau[1], j, i, dj);
        drv[3](ijk) = grad<1>(tau[2], j, i, dj);

        if ((stress_diff_t)ct_params_t::stress_diff == pade)
        {
          for (int d = 0; d < 4; ++d)
            pade_deriv(d);
        }

        switch ((sgs_scheme_t)ct_params_t::sgs_scheme)
        {
          case dns:
          {
            this->state(ix::vip_i)(ijk) += 2 * eta * this->dt * (drv[0](ijk) + drv[2](ijk));
            this->state(ix::vip_j)(ijk) += 2 * eta * this->dt * (drv[1](ijk) + drv[3](ijk));
            break;
          }
          case smg:
            assert(false);
            break;
          default:
            assert(false);
        }

        parent_t::hook_ante_step();
      }

      public:

      struct rt_params_t : parent_t::rt_params_t
      {
        typename ct_params_t::real_t eta = 0;
      };

      // ctor
      mpdata_rhs_vip_prs_sgs(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) :
	parent_t(args, p),
        eta(p.eta),
        tau(args.mem->tmp[__FILE__][0]),
        drv(args.mem->tmp[__FILE__][1]),
        wrk(args.mem->tmp[__FILE__][2])
      {
        if (eta == 0)
          throw std::runtime_error("eta == 0");
      }

      static void alloc(
        typename parent_t::mem_t *mem, 
        const int &n_iters
      ) {
	parent_t::alloc(mem, n_iters);
        parent_t::alloc_tmp_sclr(mem, __FILE__, 3); // unique strain rate tensor elements
        // TODO: do not allocate unnecessary memory when not using pade differencing
        parent_t::alloc_tmp_sclr(mem, __FILE__, 4);
        parent_t::alloc_tmp_sclr(mem, __FILE__, 2);
      }

    };

    // 3D version
    template <class ct_params_t>
    class mpdata_rhs_vip_prs_sgs<
      ct_params_t,
      typename std::enable_if<ct_params_t::n_dims == 3>::type
    > : public mpdata_rhs_vip_prs<ct_params_t>
    {
      using parent_t = mpdata_rhs_vip_prs<ct_params_t>;

      protected:

      // member fields
      arrvec_t<typename parent_t::arr_t> &tau;
      arrvec_t<typename parent_t::arr_t> &drv;
      arrvec_t<typename parent_t::arr_t> &wrk;
      typename ct_params_t::real_t eta = 0;

      template <int d, class arg_t>
      inline auto pade_helper(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k
      ) return_macro(,
	  x(idxperm::pi<d>(i+1, j, k)) + 2 * x(idxperm::pi<d>(i, j, k)) + x(idxperm::pi<d>(i - 1, j, k))
      )

      // apply pade scheme to the d-th element of the drv array
      void pade_deriv(int d)
      {
        wrk[0](this->ijk) = drv[d](this->ijk);

        for (int m = 0; m < 3; ++m)
        {
          this->xchng_sclr(wrk[0], this->ijk);
          switch(d)
          {
            case 0 : case 1 : case 2 :
            {
              wrk[1](this->ijk) = pade_helper<0>(wrk[0], this->i, this->j, this->k);
              break;
            }
            case 3 : case 4 : case 5 :
            {
              wrk[1](this->ijk) = pade_helper<1>(wrk[0], this->j, this->k, this->i);
              break;
            }
            case 6 : case 7 : case 8 :
            {
              wrk[1](this->ijk) = pade_helper<2>(wrk[0], this->k, this->i, this->j);
              break;
            }
            default : assert(false);
          }
          wrk[0](this->ijk) += (drv[d](this->ijk) - 0.25 * wrk[1](this->ijk));
        }

        drv[d](this->ijk) = wrk[0](this->ijk);
      }

      void hook_ante_step()
      {
        using namespace formulae::nabla;
        using ix = typename ct_params_t::ix;

        // to avoid writing this->
        const auto &i = this->i, &j = this->j, &k = this->k;
        const auto &di = this->di, &dj = this->dj, &dk = this->dk;
        const auto &ijk = this->ijk;

        // TODO: get rid of superfluous barriers
        this->xchng_sclr(this->state(ix::vip_i), ijk);
        this->xchng_sclr(this->state(ix::vip_j), ijk);
        this->xchng_sclr(this->state(ix::vip_k), ijk);

        // velocity gradient
        drv[0](ijk) = grad<0>(this->state(ix::vip_i), i, j, k, di);
        drv[1](ijk) = grad<0>(this->state(ix::vip_j), i, j, k, di);
        drv[2](ijk) = grad<0>(this->state(ix::vip_k), i, j, k, di);

        drv[3](ijk) = grad<1>(this->state(ix::vip_i), j, k, i, dj);
        drv[4](ijk) = grad<1>(this->state(ix::vip_j), j, k, i, dj);
        drv[5](ijk) = grad<1>(this->state(ix::vip_k), j, k, i, dj);

        drv[6](ijk) = grad<2>(this->state(ix::vip_i), k, i, j, dk);
        drv[7](ijk) = grad<2>(this->state(ix::vip_j), k, i, j, dk);
        drv[8](ijk) = grad<2>(this->state(ix::vip_k), k, i, j, dk);

        if ((stress_diff_t)ct_params_t::stress_diff == pade)
        {
          for (int d = 0; d < 9; ++d)
            pade_deriv(d);
        }

        // unique strain rate tensor elements
        tau[0](ijk) = drv[0](ijk);

        tau[1](ijk) = 0.5 * (drv[3](ijk) + drv[1](ijk));

        tau[2](ijk) = 0.5 * (drv[6](ijk) + drv[2](ijk));

        tau[3](ijk) = drv[4](ijk);

        tau[4](ijk) = 0.5 * (drv[5](ijk) + drv[7](ijk));

        tau[5](ijk) = drv[8](ijk);

        // TODO: get rid of superfluous barriers
        for (auto& t : tau)
          this->xchng_sclr(t, ijk);

        // elements of strain rate tensor divergence
        drv[0](ijk) = grad<0>(tau[0], i, j, k, di);
        drv[1](ijk) = grad<0>(tau[1], i, j, k, di);
        drv[2](ijk) = grad<0>(tau[2], i, j, k, di);

        drv[3](ijk) = grad<1>(tau[1], j, k, i, dj);
        drv[4](ijk) = grad<1>(tau[3], j, k, i, dj);
        drv[5](ijk) = grad<1>(tau[4], j, k, i, dj);

        drv[6](ijk) = grad<2>(tau[2], k, i, j, dk);
        drv[7](ijk) = grad<2>(tau[4], k, i, j, dk);
        drv[8](ijk) = grad<2>(tau[5], k, i, j, dk);

        if ((stress_diff_t)ct_params_t::stress_diff == pade)
        {
          for (int d = 0; d < 9; ++d)
            pade_deriv(d);
        }

        switch ((sgs_scheme_t)ct_params_t::sgs_scheme)
        {
          case dns:
          {
            this->state(ix::vip_i)(ijk) += 2 * eta * this->dt * (drv[0](ijk) + drv[3](ijk) + drv[6](ijk));
            this->state(ix::vip_j)(ijk) += 2 * eta * this->dt * (drv[1](ijk) + drv[4](ijk) + drv[7](ijk));
            this->state(ix::vip_k)(ijk) += 2 * eta * this->dt * (drv[2](ijk) + drv[5](ijk) + drv[8](ijk));
            break;
          }
          case smg:
            assert(false);
            break;
          default:
            assert(false);
        }

        parent_t::hook_ante_step();
      }

      public:

      struct rt_params_t : parent_t::rt_params_t
      {
        typename ct_params_t::real_t eta = 0;
      };

      // ctor
      mpdata_rhs_vip_prs_sgs(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) :
	parent_t(args, p),
        eta(p.eta),
        tau(args.mem->tmp[__FILE__][0]),
        drv(args.mem->tmp[__FILE__][1]),
        wrk(args.mem->tmp[__FILE__][2])
      {
        if (eta == 0)
          throw std::runtime_error("eta == 0");
      }

      static void alloc(
        typename parent_t::mem_t *mem,
        const int &n_iters
      ) {
	parent_t::alloc(mem, n_iters);
        parent_t::alloc_tmp_sclr(mem, __FILE__, 6); // unique strain rate tensor elements
        // TODO: do not allocate unnecessary memory when not using pade differencing
        parent_t::alloc_tmp_sclr(mem, __FILE__, 9);
        parent_t::alloc_tmp_sclr(mem, __FILE__, 2);
      }
    };
  } // namespace solvers
} // namespace libmpdataxx
