/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>

template <class ct_params_t>
class convergence_1d_spacetime : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  real_t advector(real_t t, real_t x)
  {
    return cos(t) * exp(cos(x)) / (2 + sin(x) * sin(t));
  }

  bool calc_gc() final
  {
    using libmpdataxx::opts::isset;
    using libmpdataxx::opts::div_3rd;
    using namespace libmpdataxx::arakawa_c;

    auto t = this->time;
    auto dt = this->dt;

    for (int i = this->i.first()-1; i <= this->i.last(); ++i) // starting at i.first()-1, because MPI requires that vector to the left of the domain is calculated by thread rank 0
    {
      auto x = (i+0.5) * this->di;
      this->mem->GC[0](i+h) = dt / this->di * advector(t + 0.5 * dt, x);

      if (isset(ct_params_t::opts, div_3rd))
      {
        this->mem->ndt_GC[0](i+h) = dt / this->di * (advector(t, x) - advector(t - dt, x));
        this->mem->ndtt_GC[0](i+h) = dt / this->di * (advector(t, x) + advector(t - 2 * dt, x) - 2 * advector(t - dt, x));
      }
    }

    this->xchng_vctr_alng(this->mem->GC);
    if (isset(ct_params_t::opts, div_3rd))
    {
      this->xchng_vctr_alng(this->mem->ndt_GC);
      this->xchng_vctr_alng(this->mem->ndtt_GC);
    }

    return true;
  }
};
