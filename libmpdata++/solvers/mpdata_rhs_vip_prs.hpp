/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_2d_gcrk.hpp> 
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_3d_gcrk.hpp> 
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_2d_mr.hpp> 
#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_2d_pc.hpp> 

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

    // the mpdata class
    template<typename ct_params_t, class enableif = void> 
    class mpdata_rhs_vip_prs
    {
      static_assert(!std::is_void<enableif>::value, "please specify pressure scheme type !");
    };

    // minimal residual 2D
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)mr && ct_params_t::n_dims == 2>::type
    > : public detail::mpdata_rhs_vip_prs_2d_mr<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_2d_mr<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors
    };

    // conjugate residual 2D
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)cr && ct_params_t::n_dims == 2>::type
    > : public detail::mpdata_rhs_vip_prs_2d_gcrk<ct_params_t, 1>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_2d_gcrk<ct_params_t, 1>; 
      using parent_t::parent_t; // inheriting constructors
    };
    
    // generalized conjugate residual 2D
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)gcrk && ct_params_t::n_dims == 2>::type
    > : public detail::mpdata_rhs_vip_prs_2d_gcrk<ct_params_t, ct_params_t::prs_k_iters>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_2d_gcrk<ct_params_t, ct_params_t::prs_k_iters>; 
      using parent_t::parent_t; // inheriting constructors
    };
    
    // conjugate residual 3D
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)cr && ct_params_t::n_dims == 3>::type
    > : public detail::mpdata_rhs_vip_prs_3d_gcrk<ct_params_t, 1>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_3d_gcrk<ct_params_t, 1>;
      using parent_t::parent_t; // inheriting constructors
    };

    // generalized conjugate residual 3D
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)gcrk && ct_params_t::n_dims == 3>::type
    > : public detail::mpdata_rhs_vip_prs_3d_gcrk<ct_params_t, ct_params_t::prs_k_iters>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_3d_gcrk<ct_params_t, ct_params_t::prs_k_iters>; 
      using parent_t::parent_t; // inheriting constructors
    };

    // preconditioned 2D
    template<typename ct_params_t>
    class mpdata_rhs_vip_prs<
      ct_params_t,
      typename std::enable_if<(int)ct_params_t::prs_scheme == (int)pc && ct_params_t::n_dims == 2>::type
    > : public detail::mpdata_rhs_vip_prs_2d_pc<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_prs_2d_pc<ct_params_t>; 
      using parent_t::parent_t; // inheriting constructors
    };
  }; // namespace solvers
}; // namescpae libmpdataxx
