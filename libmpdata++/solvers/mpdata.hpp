/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_osc_1d.hpp>
#include <libmpdata++/solvers/detail/mpdata_osc_2d.hpp>
#include <libmpdata++/solvers/detail/mpdata_osc_3d.hpp>

#include <libmpdata++/solvers/detail/mpdata_fct_1d.hpp>
#include <libmpdata++/solvers/detail/mpdata_fct_2d.hpp>
#include <libmpdata++/solvers/detail/mpdata_fct_3d.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    struct mpdata_family_tag {};

    // the mpdata class
    template<typename ct_params_t, int minhalo = 0, class enableif = void>
    class mpdata
    {};

    // oscillatory version
    template<typename ct_params_t, int minhalo>
    class mpdata<
      ct_params_t, minhalo,
      typename std::enable_if<!opts::isset(ct_params_t::opts, opts::fct)>::type
    > : public detail::mpdata_osc<ct_params_t, minhalo>
    {
      using parent_t = detail::mpdata_osc<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors

      protected:
      using solver_family = mpdata_family_tag;
    };

    // non-oscillatory version
    template<typename ct_params_t, int minhalo>
    class mpdata<
      ct_params_t, minhalo,
      typename std::enable_if<opts::isset(ct_params_t::opts, opts::fct)>::type
    > : public detail::mpdata_fct<ct_params_t, minhalo>
    {
      using parent_t = detail::mpdata_fct<ct_params_t, minhalo>;
      using parent_t::parent_t; // inheriting constructors

      protected:
      using solver_family = mpdata_family_tag;
    };
  } // namespace solvers
} // namescpae libmpdataxx
