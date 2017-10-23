/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>

template <class ct_params_t>
class convergence_2d_spacetime : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  real_t advector(real_t t, real_t x, real_t y, real_t comp)
  {
    return cos(t) * exp(cos(x) + cos(y)) / (2 + sin(comp) * sin(t));
  }

  bool calc_gc() final
  {
    using libmpdataxx::opts::isset;
    using libmpdataxx::opts::div_3rd;
    using namespace libmpdataxx::arakawa_c;

    auto t = this->time;
    auto dt = this->dt;

    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        auto x = (i+0.5) * this->di;
        auto y = j * this->dj;
        this->mem->GC[0](i+h, j) = dt / this->di * advector(t + 0.5 * dt, x, y, x);

        if (isset(ct_params_t::opts, div_3rd))
        {
          this->mem->ndt_GC[0](i+h , j) = dt / this->di *
                                          (advector(t, x, y, x) - advector(t - dt, x, y, x));
          this->mem->ndtt_GC[0](i+h, j) = dt / this->di *
                                          (advector(t, x, y, x) + advector(t - 2 * dt, x, y, x) - 2 * advector(t - dt, x, y, x));
        }
      }
    }
    
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        auto x = i * this->di;
        auto y = (j+0.5) * this->dj;
        this->mem->GC[1](i, j+h) = dt / this->dj * advector(t + 0.5 * dt, x, y, y);

        if (isset(ct_params_t::opts, div_3rd))
        {
          this->mem->ndt_GC[1](i, j+h ) = dt / this->dj *
                                          (advector(t, x, y, y) - advector(t - dt, x, y, y));
          this->mem->ndtt_GC[1](i, j+h) = dt / this->dj *
                                          (advector(t, x, y, y) + advector(t - 2 * dt, x, y, y) - 2 * advector(t - dt, x, y, y));
        }
      }
    }

    auto ex = this->halo - 1;
    this->xchng_vctr_alng(this->mem->GC);
    this->xchng_vctr_nrml(this->mem->GC, this->ijk, ex);
    
    if (isset(ct_params_t::opts, div_3rd))
    {
      this->xchng_vctr_alng(this->mem->ndt_GC);
      this->xchng_vctr_nrml(this->mem->ndt_GC, this->ijk, ex);
      
      this->xchng_vctr_alng(this->mem->ndtt_GC);
      this->xchng_vctr_nrml(this->mem->ndtt_GC, this->ijk, ex);
    }
    return true;
  }
};
