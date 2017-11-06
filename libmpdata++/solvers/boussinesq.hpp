/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/boussinesq_expl.hpp>
#include <libmpdata++/solvers/detail/boussinesq_impl.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    struct mpdata_boussinesq_family_tag {};
    struct mpdata_boussinesq_sgs_family_tag {};

    // the boussinesq class
    template<typename ct_params_t, class enableif = void>
    class boussinesq
    {};

    template<typename ct_params_t>
    class boussinesq<
      ct_params_t,
      typename std::enable_if<!ct_params_t::impl_tht>::type
    > : public detail::boussinesq_expl<ct_params_t>
    {
      using parent_t = detail::boussinesq_expl<ct_params_t>;
      using parent_t::parent_t; // inheriting constructors

      protected:
      using solver_family = typename std::conditional<static_cast<sgs_scheme_t>(ct_params_t::sgs_scheme) == iles,
                                                                  mpdata_boussinesq_family_tag,
                                                                  mpdata_boussinesq_sgs_family_tag>::type;
    };
    
    template<typename ct_params_t>
    class boussinesq<
      ct_params_t,
      typename std::enable_if<ct_params_t::impl_tht>::type
    > : public detail::boussinesq_impl<ct_params_t>
    {
      using parent_t = detail::boussinesq_impl<ct_params_t>;
      using parent_t::parent_t; // inheriting constructors
      
      protected:
      using solver_family = typename std::conditional<static_cast<sgs_scheme_t>(ct_params_t::sgs_scheme) == iles,
                                                      mpdata_boussinesq_family_tag,
                                                      mpdata_boussinesq_sgs_family_tag>::type;
    };
  } // namespace solvers
} // namescpae libmpdataxx
