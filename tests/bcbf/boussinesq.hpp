/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>

template <class ct_params_t>
class boussinesq : public libmpdataxx::solvers::mpdata_rhs_vip_prs<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata_rhs_vip_prs<ct_params_t>;
  using ix = typename ct_params_t::ix;

  public:
  using real_t = typename ct_params_t::real_t;

  private:
  // member fields
  real_t g, Tht_ref;

//<listing-1>
  // explicit forcings 
  void update_rhs(
    libmpdataxx::arrvec_t<
      typename parent_t::arr_t
    > &rhs, 
    const real_t &dt, 
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at); 

    const auto Tht  = this->psi_n(ix::tht); 
    const auto &ijk = this->ijk;

    rhs.at(ix::w)(ijk) += 
      g * (Tht(ijk) - Tht_ref) / Tht_ref; 
  }
//</listing-1>
//    rhs.at(ix::w)(ijk) += g /*/ 300*/ * (Tht(ijk) - Tht_ref) / Tht_ref; 

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t g = 9.81, Tht_ref = 0; 
  };

  // ctor
  boussinesq( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p),
    g(p.g),
    Tht_ref(p.Tht_ref)
  {
    assert(Tht_ref != 0);
  }
};
