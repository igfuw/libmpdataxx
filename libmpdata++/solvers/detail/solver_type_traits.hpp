#pragma once

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      // helpers that detect if the solver uses fractal reconstruction
      template< class, class = void >
      struct slvr_with_frac_recn : std::false_type { };

      template< class ct_params_t >
      struct slvr_with_frac_recn<
        ct_params_t, 
        typename std::enable_if_t<(int)ct_params_t::fractal_recon != (int)0 > > : std::true_type { };
        //std::enable_if_t< std::is_base_of_v< libmpdataxx::solvers::mpdata_rhs_vip_prs_sgs_fra_family_tag, solver_t > > > : std::true_type { };

      template <class ct_params_t>
      inline constexpr bool slvr_with_frac_recn_v = slvr_with_frac_recn<ct_params_t>::value;
    }
  }
}
