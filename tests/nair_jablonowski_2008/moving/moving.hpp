/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>
#include "../common/transforms.hpp"

template <class ct_params_t>
class moving : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  real_t a, u0, v0, x0, y0, xc, yc;

  protected:

  xpf_t xpf;
  ixpf_t ixpf;
  ypf_t ypf;
  iypf_t iypf;

  void update_vortex_centre(real_t t)
  {
    xpf.x0 = pi;
    xpf.y0 = pi / 2 - a;
    ypf.x0 = pi;
    ypf.y0 = pi / 2 - a;

    ixpf.x0 = pi;
    ixpf.y0 = pi / 2 - a;
    iypf.x0 = pi;
    iypf.y0 = pi / 2 - a;

    real_t x = xpf(x0, y0);
    real_t y = ypf(x0, y0);

    x += u0 * t;

    xc = ixpf(x, y);
    yc = iypf(x, y);
    
    xpf.x0 = xc;
    xpf.y0 = yc;
    ypf.x0 = xc;
    ypf.y0 = yc;
  }

  real_t vel_x(real_t x, real_t y)
  {
    real_t r = 3 * cos(ypf(x, y));
    real_t omg = r != 0 ? v0 * 3 * sqrt(2.) / (2 * r) * tanh(r) / pow(cosh(r), 2) : 0;
    return (   omg * (sin(yc) - cos(yc) * cos(x - xc) * tan(y))
             + u0 * (cos(a) + tan(y) * sin(a) * cos(x))
           ) * this->di * this->dj * cos(y) * this->dt / this->di;
  }
  
  real_t vel_y(real_t x, real_t y)
  {
    real_t r = 3 * cos(ypf(x, y));
    real_t omg = r != 0 ? v0 * 3 * sqrt(2.) / (2 * r) * tanh(r) / pow(cosh(r), 2) : 0;
    return (   omg * cos(yc) * sin(x - xc)
             - u0 * sin(x) * sin(a)
           ) * this->di * this->dj * cos(y) * this->dt/ this->dj;
  }

  bool calc_gc() final
  {
    // commented out code will be used when my new mpdata version lands on master !

    //using libmpdataxx::opts::isset;
    //using libmpdataxx::opts::amz;

    using namespace libmpdataxx::arakawa_c;
    auto t = this->time;
    //if (isset(ct_params_t::opts, amz))
    //{
    //  update_vortex_centre(t - 0.5 * this->dt);
    //  for (int i = this->i.first(); i <= this->i.last(); ++i)
    //  {
    //    for (int j = this->j.first(); j <= this->j.last(); ++j)
    //    {
    //      auto x = (i+0.5) * this->di;
    //      auto y = (j+0.5) * this->dj - pi / 2;
    //      this->mem->dGC_dt[0](i+h, j) = -vel_x(x, y);
    //      this->mem->dGC_dtt[0](i+h, j) = vel_x(x, y);

    //      x = i * this->di;
    //      y = (j+1) * this->dj - pi / 2;
    //      this->mem->dGC_dt[1](i, j+h) = -vel_y(x, y);
    //      this->mem->dGC_dtt[1](i, j+h) = vel_y(x, y);
    //    }
    //  }
    //  update_vortex_centre(t);
    //  for (int i = this->i.first(); i <= this->i.last(); ++i)
    //  {
    //    for (int j = this->j.first(); j <= this->j.last(); ++j)
    //    {
    //      auto x = (i+0.5) * this->di;
    //      auto y = (j+0.5) * this->dj - pi / 2;
    //      this->mem->dGC_dtt[0](i+h, j) -= 2 * vel_x(x, y);

    //      x = i * this->di;
    //      y = (j+1) * this->dj - pi / 2;
    //      this->mem->dGC_dtt[1](i, j+h) -= 2 * vel_y(x, y);
    //    }
    //  }
    //}

    update_vortex_centre(t + 0.5 * this->dt);
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        auto x = (i+0.5) * this->di;
        auto y = (j+0.5) * this->dj - pi / 2;
        this->mem->GC[0](i+h, j) = vel_x(x, y);
        //if (isset(ct_params_t::opts, amz))
        //{
        //  this->mem->dGC_dt[0](i+h, j) += vel_x(x, y);
        //  this->mem->dGC_dtt[0](i+h, j) += vel_x(x, y);
        //}

        x = i * this->di;
        y = (j+1) * this->dj - pi / 2;
        this->mem->GC[1](i, j+h) = vel_y(x, y);
        //if (isset(ct_params_t::opts, amz))
        //{
        //  this->mem->dGC_dt[1](i, j+h) += vel_y(x, y);
        //  this->mem->dGC_dtt[1](i, j+h) += vel_y(x, y);
        //}
      }
    }

    auto ex = this->halo - 1;
    this->xchng_vctr_alng(this->mem->GC);
    this->xchng_vctr_nrml(this->mem->GC, this->i^ex, this->j^ex);
    
    //if (isset(ct_params_t::opts, amz))
    //{
    //  this->xchng_vctr_alng(this->mem->dGC_dt);
    //  this->xchng_vctr_nrml(this->mem->dGC_dt, this->i^ex, this->j^ex);
    //  
    //  this->xchng_vctr_alng(this->mem->dGC_dtt);
    //  this->xchng_vctr_nrml(this->mem->dGC_dtt, this->i^ex, this->j^ex);
    //}

    return true;
  }

  public:
  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t x0, y0, a, u0, v0;
  };

  // ctor
  moving( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p),
    a(p.a),
    u0(p.u0),
    v0(p.v0),
    x0(p.x0),
    y0(p.y0),
    xpf{.x0 = x0, .y0 = y0},
    ypf{.x0 = x0, .y0 = y0}
  {}
};
