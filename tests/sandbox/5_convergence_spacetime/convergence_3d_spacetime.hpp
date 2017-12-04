/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>

template <class ct_params_t>
class convergence_3d_spacetime : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  real_t advector(real_t t, real_t x, real_t y, real_t z, real_t comp)
  {
    return cos(t) * exp(cos(x) + cos(y) + cos(z)) / (2 + sin(comp) * sin(t));
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
        for (int k = this->k.first(); k <= this->k.last(); ++k)
        {
          auto x = (i+0.5) * this->di;
          auto y = j * this->dj;
          auto z = k * this->dk;

          this->mem->GC[0](i+h, j, k) = dt / this->di * advector(t + 0.5 * dt, x, y, z, x);

          if (isset(ct_params_t::opts, div_3rd))
          {
            this->mem->ndt_GC[0](i+h, j, k ) = dt / this->di *
                                            (advector(t, x, y, z, x) - advector(t - dt, x, y, z, x));
            this->mem->ndtt_GC[0](i+h, j, k) = dt / this->di * (
                                               advector(t, x, y, z, x)
                                             + advector(t - 2 * dt, x, y, z, x)
                                             - 2 * advector(t - dt, x, y, z, x)
                                             );
          }
        }
      }
    }
    
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        for (int k = this->k.first(); k <= this->k.last(); ++k)
        {
          auto x = i * this->di;
          auto y = (j+0.5) * this->dj;
          auto z = k * this->dk;

          this->mem->GC[1](i, j+h, k) = dt / this->dj * advector(t + 0.5 * dt, x, y, z, y);

          if (isset(ct_params_t::opts, div_3rd))
          {
            this->mem->ndt_GC[1](i, j+h, k ) = dt / this->dj *
                                            (advector(t, x, y, z, y) - advector(t - dt, x, y, z, y));
            this->mem->ndtt_GC[1](i, j+h, k) = dt / this->dj * (
                                               advector(t, x, y, z, y)
                                             + advector(t - 2 * dt, x, y, z, y)
                                             - 2 * advector(t - dt, x, y, z, y)
                                             );
          }
        }
      }
    }
    
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        for (int k = this->k.first(); k <= this->k.last(); ++k)
        {
          auto x = i * this->di;
          auto y = j * this->dj;
          auto z = (k+0.5) * this->dk;

          this->mem->GC[2](i, j, k+h) = dt / this->dk * advector(t + 0.5 * dt, x, y, z, z);

          if (isset(ct_params_t::opts, div_3rd))
          {
            this->mem->ndt_GC[2](i, j, k+h ) = dt / this->dk *
                                            (advector(t, x, y, z, z) - advector(t - dt, x, y, z, z));
            this->mem->ndtt_GC[2](i, j, k+h) = dt / this->dk * (
                                               advector(t, x, y, z, z)
                                             + advector(t - 2 * dt, x, y, z, z)
                                             - 2 * advector(t - dt, x, y, z, z)
                                             );
          }
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
