/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template <class ct_params_t>
      class boussinesq_common : public libmpdataxx::solvers::mpdata_rhs_vip_prs<ct_params_t>
      {
        using parent_t = libmpdataxx::solvers::mpdata_rhs_vip_prs<ct_params_t>;

        public:
        using real_t = typename ct_params_t::real_t;

        protected:
        // member fields
        real_t g, Tht_ref;
        typename parent_t::arr_t &tht_e, &tht_abs, &hflux;

        public:
        struct rt_params_t : parent_t::rt_params_t 
        { 
          real_t g = 9.81, Tht_ref = 0; 
        };

        // ctor
        boussinesq_common( 
          typename parent_t::ctor_args_t args, 
          const rt_params_t &p
        ) :
          parent_t(args, p),
          g(p.g),
          Tht_ref(p.Tht_ref),
          tht_e(args.mem->tmp[__FILE__][0][0]),
          tht_abs(args.mem->tmp[__FILE__][1][0]),
          hflux(args.mem->tmp[__FILE__][2][0])
        {
          assert(Tht_ref != 0);
          hflux(this->ijk) = 0;
          tht_abs(this->ijk) = 0;
        }

        static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
        {
          parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1, "tht_e");
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1, "tht_abs");
          parent_t::alloc_tmp_sclr(mem, __FILE__, 1, "hflux");
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
