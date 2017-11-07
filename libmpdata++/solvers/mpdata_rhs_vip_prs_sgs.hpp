/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_dns.hpp>
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_smg.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    enum sgs_scheme_t
    {
      iles,
      dns,
      smg
    };

    const std::map<sgs_scheme_t, std::string> sgs2string = {
      {iles, "iles"},
      {dns , "dns" },
      {smg , "smg" }
    };

    struct mpdata_rhs_vip_prs_sgs_family_tag {};
    struct mpdata_rhs_vip_prs_sgs_dns_family_tag {};
    struct mpdata_rhs_vip_prs_sgs_smg_family_tag {};

    template<typename ct_params_t, class enableif = void> 
    class mpdata_rhs_vip_prs_sgs;

    template<typename ct_params_t>
    class mpdata_rhs_vip_prs_sgs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::sgs_scheme == (int)iles>::type
    > : public mpdata_rhs_vip_prs<ct_params_t>
    {
      using parent_t = mpdata_rhs_vip_prs<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors
    };
    
    template <class ct_params_t>
    class mpdata_rhs_vip_prs_sgs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::sgs_scheme == (int)dns>::type
    > : public detail::mpdata_rhs_vip_prs_sgs_dns<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_sgs_dns<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors

      protected:
      using solver_family = mpdata_rhs_vip_prs_sgs_dns_family_tag;
    };
    
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs_sgs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::sgs_scheme == (int)smg>::type
    > : public detail::mpdata_rhs_vip_prs_sgs_smg<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_sgs_smg<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors
      
      protected:
      using solver_family = mpdata_rhs_vip_prs_sgs_smg_family_tag;
    };
  } // namespace solvers
} // namespace libmpdataxx
