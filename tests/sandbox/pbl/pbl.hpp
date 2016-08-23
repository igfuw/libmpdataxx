/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/boussinesq.hpp>
#include <cmath>

template <class ct_params_t>
class pbl : public libmpdataxx::solvers::boussinesq<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::boussinesq<ct_params_t>;
  using ix = typename ct_params_t::ix;

  public:
  using real_t = typename ct_params_t::real_t;

  private:
  real_t hscale, cdrag;

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();

    for (int k = this->k.first(); k <= this->k.last(); ++k)
    {
      this->vip_rhs[0](this->i, this->j, k) += - cdrag / hscale * this->dt * sqrt(
                                                pow2(this->state(ix::vip_i)(this->i, this->j, 0))
                                              + pow2(this->state(ix::vip_j)(this->i, this->j, 0))
                                              ) * this->state(ix::vip_i)(this->i, this->j, 0)
                                                * exp(-this->dj * k / hscale);
      
      this->vip_rhs[1](this->i, this->j, k) += - cdrag / hscale * this->dt * sqrt(
                                                pow2(this->state(ix::vip_i)(this->i, this->j, 0))
                                              + pow2(this->state(ix::vip_j)(this->i, this->j, 0))
                                              ) * this->state(ix::vip_j)(this->i, this->j, 0)
                                                * exp(-this->dj * k / hscale);
    }
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t hscale = 1, cdrag = 0; 
  };

  // ctor
  pbl( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p),
    hscale(p.hscale),
    cdrag(p.cdrag)
  {}
};
