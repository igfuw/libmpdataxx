/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_gcrk.hpp> 
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
      gcrk, // generalized conjugate residual (restarted after k steps)
      pc  // preconditioned
    };

    const std::map<prs_scheme_t, std::string> prs2string = {
      {mr, "mr"},
      {cr, "cr"},
      {gcrk, "gcrk"},
      {pc, "pc"}
    };

    struct mpdata_rhs_vip_prs_family_tag {};

    // the mpdata class
    template<typename ct_params_t, int minhalo = 0, class enableif = void> 
    class mpdata_rhs_vip_prs
    {
      static_assert(!std::is_void<enableif>::value, "please specify pressure scheme type !");
    };

    // minimal residual
    template<typename ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs<
      ct_params_t, minhalo,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)mr>::type
    > : public detail::mpdata_rhs_vip_prs_mr<ct_params_t, minhalo>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_mr<ct_params_t, minhalo>; 
      using parent_t::parent_t; // inheriting constructors

      protected:
      using solver_family = mpdata_rhs_vip_prs_family_tag;
    };

    // conjugate residual
    template<typename ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs<
      ct_params_t, minhalo,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)cr>::type
    > : public detail::mpdata_rhs_vip_prs_gcrk<ct_params_t, 1, minhalo>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_gcrk<ct_params_t, 1, minhalo>; 
      using parent_t::parent_t; // inheriting constructors
      
      protected:
      using solver_family = mpdata_rhs_vip_prs_family_tag;
    };
    
    // generalized conjugate residual
    template<typename ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs<
      ct_params_t, minhalo,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)gcrk>::type
    > : public detail::mpdata_rhs_vip_prs_gcrk<ct_params_t, ct_params_t::prs_k_iters, minhalo>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_gcrk<ct_params_t, ct_params_t::prs_k_iters, minhalo>; 
      using parent_t::parent_t; // inheriting constructors
      
      protected:
      using solver_family = mpdata_rhs_vip_prs_family_tag;
    };

    // preconditioned
    template<typename ct_params_t, int minhalo>
    class mpdata_rhs_vip_prs<
      ct_params_t, minhalo,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)pc>::type
    > : public detail::mpdata_rhs_vip_prs_pc<ct_params_t, minhalo>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_pc<ct_params_t, minhalo>; 
      using parent_t::parent_t; // inheriting constructors
      
      protected:
      using solver_family = mpdata_rhs_vip_prs_family_tag;
    };
  } // namespace solvers
} // namescpae libmpdataxx
