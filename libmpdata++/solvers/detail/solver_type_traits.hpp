#pragma once

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      // helper function that detect if the solver uses fractal reconstruction
      template< class, class = void >
      struct slvr_with_frac_recn : std::false_type { };

      template< class solver_t >
      struct slvr_with_frac_recn<solver_t, std::void_t<typename solver_t::mpdata_rhs_vip_prs_sgs_fra_family_tag>> : std::true_type { };
    }
  }
}
