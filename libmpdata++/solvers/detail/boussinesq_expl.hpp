/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include <libmpdata++/solvers/detail/boussinesq_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class boussinesq_expl : public boussinesq_common<ct_params_t>
      {
        using parent_t = boussinesq_common<ct_params_t>;

        public:
        using real_t = typename ct_params_t::real_t;

        protected:
        // member fields
        using ix = typename ct_params_t::ix;
        bool buoy_filter;
        typename parent_t::arr_t &tmp1, &tmp2;
        
        template <int nd = ct_params_t::n_dims> 
        void filter(typename std::enable_if<nd == 2>::type* = 0)
        {
          const auto &i(this->i), &j(this->j);
          this->xchng_sclr(tmp1, i, j);
          tmp2(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
        }

        template <int nd = ct_params_t::n_dims> 
        void filter(typename std::enable_if<nd == 3>::type* = 0)
        {
          const auto &i(this->i), &j(this->j), &k(this->k);
          this->xchng_sclr(tmp1, i, j, k);
          tmp2(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k + 1));
        }
        
        // helpers for buoyancy forces
        template<class ijk_t>
        inline auto buoy_at_0(const ijk_t &ijk)
        {
          return return_helper<rng_t>(
            this->g * (this->state(ix::tht)(ijk) - this->tht_e(ijk)) / this->Tht_ref
          );
        }
        
        template<class ijk_t>
        inline auto buoy_at_1(const ijk_t &ijk)
        {
          return return_helper<rng_t>(
            this->g * (
                (  this->state(ix::tht)(ijk)
                 + 0.5 * this->dt * this->hflux(ijk) + 0.5 * this->dt * this->tht_abs(ijk) * this->tht_e(ijk))
                / (1 + 0.5 * this->dt * this->tht_abs(ijk))
                - this->tht_e(ijk)
              ) / this->Tht_ref
          );
        }

        // explicit forcings 
        void update_rhs(
          libmpdataxx::arrvec_t<
            typename parent_t::arr_t
          > &rhs, 
          const real_t &dt, 
          const int &at 
        ) {
          parent_t::update_rhs(rhs, dt, at); 

          const auto &tht = this->state(ix::tht); 
          const auto &ijk = this->ijk;

          auto ix_w = this->vip_ixs[ct_params_t::n_dims - 1];

          switch (at)
          {
            case (0):
            {
              rhs.at(ix::tht)(ijk) += this->hflux(ijk) - this->tht_abs(ijk) * (tht(ijk) - this->tht_e(ijk));
              
              if (!buoy_filter)
              {
                rhs.at(ix_w)(ijk) += buoy_at_0(ijk);
              }
              else
              {
                tmp1(ijk) = buoy_at_0(ijk);
                filter();
                rhs.at(ix_w)(ijk) += (tmp2)(ijk);
              }
              break;
            }
            case (1):
            {
              rhs.at(ix::tht)(ijk) += this->hflux(ijk) - this->tht_abs(ijk) * (
                (tht(ijk) + 0.5 * this->dt * this->hflux(ijk) + 0.5 * this->dt * this->tht_abs(ijk) * this->tht_e(ijk))
                / (1 + 0.5 * this->dt * this->tht_abs(ijk))
                - this->tht_e(ijk)
              );

              if (!buoy_filter)
              {
                rhs.at(ix_w)(ijk) += buoy_at_1(ijk);
              }
              else
              {
                tmp1(ijk) = buoy_at_1(ijk);
                filter();
                rhs.at(ix_w)(ijk) += (tmp2)(ijk);
              }
            }
          }
        }
        
        public:
        struct rt_params_t : parent_t::rt_params_t 
        { 
          bool buoy_filter = false; 
        };

        // ctor
        boussinesq_expl( 
          typename parent_t::ctor_args_t args, 
          const rt_params_t &p
        ) :
          parent_t(args, p),
          buoy_filter(p.buoy_filter),
          tmp1(args.mem->tmp[__FILE__][0][0]),
          tmp2(args.mem->tmp[__FILE__][0][1])
        {}

        static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
        {
          parent_t::alloc(mem, n_iters);
          // TODO: do not allocate if not filtering
          parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // tmp1, tmp2
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
