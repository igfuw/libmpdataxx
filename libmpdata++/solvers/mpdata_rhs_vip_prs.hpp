/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_cr.hpp> 
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_mr.hpp> 
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_pc.hpp> 

namespace libmpdataxx
{
  namespace solvers
  {
    enum prs_scheme_t 
    {   
      mr, // minimal residual
      cr, // conjugate residual
      pc  // preconditioned
    };  

    // the mpdata class
    template<typename ct_params_t, class enableif = void> 
    class mpdata_rhs_vip_prs
    {};

    // minimal residual 
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<ct_params_t::prs_scheme == mr>::type
    > : public detail::mpdata_rhs_vip_prs_mr<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_mr<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors
    };

    // conjugate residual 
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<ct_params_t::prs_scheme == cr>::type
    > : public detail::mpdata_rhs_vip_prs_cr<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_cr<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors
    };

    // preconditioned
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<ct_params_t::prs_scheme == pc>::type
    > : public detail::mpdata_rhs_vip_prs_pc<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_pc<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors
    };
  }; // namespace solvers
}; // namescpae libmpdataxx
