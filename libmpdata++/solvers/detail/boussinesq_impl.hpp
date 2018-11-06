/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include <libmpdata++/solvers/detail/boussinesq_sgs_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class boussinesq_impl : public
                              std::conditional<static_cast<sgs_scheme_t>(ct_params_t::sgs_scheme) == iles,
                                               boussinesq_common<ct_params_t>,
                                               boussinesq_sgs_common<ct_params_t>>::type
      {
        using parent_t = typename std::conditional<static_cast<sgs_scheme_t>(ct_params_t::sgs_scheme) == iles,
                                                   boussinesq_common<ct_params_t>,
                                                   boussinesq_sgs_common<ct_params_t>>::type;
        public:
        using real_t = typename ct_params_t::real_t;

        protected:
        // member fields
        using ix = typename ct_params_t::ix;
        typename parent_t::arr_t &dtht_e;
        
        template <int nd = ct_params_t::n_dims>
        void calc_dtht_e(typename std::enable_if<nd == 2>::type* = 0)
        {
          this->xchng_sclr(this->tht_e, this->ijk);
          this->dtht_e(this->ijk) = formulae::nabla::grad<1>(this->tht_e, this->j, this->i, this->dj);
        }
        
        template <int nd = ct_params_t::n_dims>
        void calc_dtht_e(typename std::enable_if<nd == 3>::type* = 0)
        {
          this->xchng_sclr(this->tht_e, this->ijk);
          this->dtht_e(this->ijk) = formulae::nabla::grad<2>(this->tht_e, this->k, this->i, this->j, this->dk);
        }

        void calc_full_tht(typename parent_t::arr_t &full_tht) final
        {
          full_tht(this->ijk) = this->state(ix::tht)(this->ijk) + this->tht_e(this->ijk);
        }

	void hook_ante_loop(const typename parent_t::advance_arg_t nt) 
	{   
          calc_dtht_e();
	  parent_t::hook_ante_loop(nt);
	}

        virtual void normalize_vip(const arrvec_t<typename parent_t::arr_t> &v)
        {
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl)
          {
            for (int d = 0; d < ct_params_t::n_dims - 1; ++d)
            {
              v[d](this->ijk) /= (1 + 0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk));
            }
            v[ct_params_t::n_dims - 1](this->ijk) /=
            (1 + 0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk)
               + 0.25 * this->dt * this->dt * this->g / this->Tht_ref * this->dtht_e(this->ijk)
                 / (1 + 0.5 * this->dt * this->tht_abs(this->ijk)));
          }
          else
          {
            v[ct_params_t::n_dims - 1](this->ijk) /=
            (1 + 0.25 * this->dt * this->dt * this->g / this->Tht_ref * this->dtht_e(this->ijk)
                 / (1 + 0.5 * this->dt * this->tht_abs(this->ijk)));
          }
        }
        
        void update_rhs(
          libmpdataxx::arrvec_t<
            typename parent_t::arr_t
          > &rhs, 
          const real_t &dt, 
          const int &at 
        ) {
          parent_t::update_rhs(rhs, dt, at);

          auto ix_w = this->vip_ixs[ct_params_t::n_dims - 1];

          const auto &tht = this->state(ix::tht); 
          const auto &w = this->state(ix_w);
          const auto &ijk = this->ijk;

          switch (at)
          {
            case (0):
            {
              rhs.at(ix::tht)(ijk) += -w(ijk) * this->dtht_e(ijk) + this->hflux_frc(ijk) - this->tht_abs(ijk) * tht(ijk);

              rhs.at(ix_w)(ijk) += this->g * tht(ijk) / this->Tht_ref;
              break;
            }
            case (1):
            {
              rhs.at(ix::tht)(ijk) += this->hflux_frc(ijk);

              rhs.at(ix_w)(ijk) += this->g * (tht(ijk) + 0.5 * this->dt * rhs.at(ix::tht)(ijk))
                                   / (this->Tht_ref * (1 + 0.5 * this->dt * this->tht_abs(ijk)));
              break;
            }
          }
        }
        
        void vip_rhs_impl_fnlz()
        {
          parent_t::vip_rhs_impl_fnlz();
          
          const auto &w = this->vips()[ct_params_t::n_dims - 1];
          this->state(ix::tht)(this->ijk) = ( this->state(ix::tht)(this->ijk) 
                                            - 0.5 * this->dt * w(this->ijk) * this->dtht_e(this->ijk))
                                            / (1 + 0.5 * this->dt * this->tht_abs(this->ijk));
          this->rhs.at(ix::tht)(this->ijk) += -w(this->ijk) * this->dtht_e(this->ijk)
                                              -this->tht_abs(this->ijk) * this->state(ix::tht)(this->ijk);
        }
        
        public:
        // ctor
        boussinesq_impl( 
          typename parent_t::ctor_args_t args, 
          const typename parent_t::rt_params_t &p
        ) :
          parent_t(args, p),
          dtht_e(args.mem->tmp[__FILE__][0][0])
        {}

        static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
        {
          parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1); // dtht_e
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
