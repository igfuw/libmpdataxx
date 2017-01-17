/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>

template <class ct_params_t>
class var_dt_test : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  real_t pi, u_0, tau;

  bool calc_gc() final
  {
    auto t_half = this->time + 0.5 * this->dt;
    using namespace libmpdataxx::arakawa_c;
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      this->mem->GC[0](i+h) = this->dt / this->di * u_0 * cos(pi * t_half / tau);
    }
    this->xchng_vctr_alng(this->mem->GC);
    return true;
  }
  
  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t pi, u_0, tau;
  };

  // ctor
  var_dt_test( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p),
    pi(p.pi),
    u_0(p.u_0),
    tau(p.tau)
  {}
};
