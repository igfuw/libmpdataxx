/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>

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
      typename std::enable_if<(int)ct_params_t::n_fra_rec == 0>::type
    > : public mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
    {
      using parent_t = mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors
    };

    template <class ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs_sgs_fra<
      ct_params_t, minhalo,
      typename std::enable_if<(int)ct_params_t::n_fra_rec > 0>::type
    > : public detail::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors

      protected:
      using solver_family = mpdata_rhs_vip_prs_sgs_fra_family_tag;
    };
  } // namespace solvers
} // namespace libmpdataxx
