/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/boussinesq_expl_2d.hpp>
#include <libmpdata++/solvers/detail/boussinesq_expl_3d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    // the boussinesq class
    template<typename ct_params_t, class enableif = void>
    class boussinesq
    {};

    template<typename ct_params_t>
    class boussinesq<
      ct_params_t,
      typename std::enable_if<ct_params_t::n_dims == 2>::type
    > : public detail::boussinesq_expl_2d<ct_params_t>
    {
      using parent_t = detail::boussinesq_expl_2d<ct_params_t>;
      using parent_t::parent_t; // inheriting constructors
    };
    
    template<typename ct_params_t>
    class boussinesq<
      ct_params_t,
      typename std::enable_if<ct_params_t::n_dims == 3>::type
    > : public detail::boussinesq_expl_3d<ct_params_t>
    {
      using parent_t = detail::boussinesq_expl_3d<ct_params_t>;
      using parent_t::parent_t; // inheriting constructors
    };
  } // namespace solvers
} // namescpae libmpdataxx
