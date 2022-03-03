/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_fra.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    struct mpdata_rhs_vip_prs_sgs_fra_family_tag {};

    template<typename ct_params_t, int minhalo = 0, class enableif = void>
    class mpdata_rhs_vip_prs_sgs_fra;

    template<typename ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs_sgs_fra<
      ct_params_t, minhalo,
      typename std::enable_if_t<(int)ct_params_t::fractal_recon == (int)0> // no fields with fractal reconstruction
//      typename std::enable_if_t<(int)ct_params_t::fractal_recon == (int) 0>
    //  std::enable_if_t<(int)ct_params_t::fractal_recon == 0, bool> = true
    > : public mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
    {
      using parent_t = mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors
    };

    template <class ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs_sgs_fra<
      ct_params_t, minhalo,
      typename std::enable_if_t<(int)ct_params_t::fractal_recon != (int)0>
      //typename std::enable_if_t<(int)ct_params_t::fractal_recon != (int) 0>
    //  std::enable_if_t<(int)ct_params_t::fractal_recon > 0, bool> = true
    > : public detail::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>,
        public mpdata_rhs_vip_prs_sgs_fra_family_tag        // allows checking if derived solvers use fractal reconstruction
    {
      using parent_t = detail::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors

      private:
      using solver_family = mpdata_rhs_vip_prs_sgs_fra_family_tag;
    };

//    template <class ct_params_t, int minhalo = 0>
//    class mpdata_rhs_vip_prs_sgs_fra : public detail::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>
//    {
//      using parent_t = detail::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>;
//      using parent_t::parent_t; // inheriting constructors
//
//      protected:
//      using solver_family = mpdata_rhs_vip_prs_sgs_fra_family_tag;
//    };
  } // namespace solvers
} // namespace libmpdataxx
