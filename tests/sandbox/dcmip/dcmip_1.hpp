/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>

template <class ct_params_t>
class dcmip_1 : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  real_t pi, a, g, z_top, p_top, p_0, rho_0, tau, H, omg_0, b;

  protected:

  static constexpr bool needs_dt = libmpdataxx::opts::isset(ct_params_t::opts, libmpdataxx::opts::div_3rd) ||
                                   libmpdataxx::opts::isset(ct_params_t::opts, libmpdataxx::opts::div_3rd_dt);
  bool calc_gc() final
  {
    using libmpdataxx::opts::isset;
    using libmpdataxx::opts::div_3rd;

    using namespace libmpdataxx::arakawa_c;
    auto t = this->time + 0.5 * this->dt;

    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        for (int k = this->k.first(); k <= this->k.last(); ++k)
        {
          auto x = i * this->di;
          auto x1 = x + 0.5 * this->di;
          auto xp = x - 2 * pi * t / tau;
          auto xp1 = xp + 0.5 * this->di;

          auto y = (j+0.5) * this->dj - pi / 2;
          auto y1 = y + 0.5 * this->dj;

          auto z = k * this->dk;
          auto z1 = z + 0.5 * this->dk;

          auto rho  = rho_0 * exp(-z / H);
          auto p  = p_0 * exp(-z / H);
          auto p1 = p_0 * exp(-z1 / H);
          auto s1 = 1 + exp((p_top - p_0) / (b * p_top)) - exp((p1 - p_0) / (b * p_top)) - exp((p_top - p1) / (b * p_top));

          auto udv = omg_0 * a / (b * p_top) * cos(xp1) * pow(cos(y), 2)
                    * (-exp((p - p_0) / (b * p_top)) + exp((p_top - p) / (b * p_top)));
                    
          auto ud = udv * cos(2 * pi * t / tau);

          this->mem->GC[0](i+h, j, k) = this->dt / this->di * 
                                        (
                                         10 * a / tau * pow(sin(xp1), 2) * sin(2 * y) * cos(pi * t / tau)
                                         + 2 * pi * a / tau * cos(y)
                                         + ud
                                        ) * rho;
          
          this->mem->GC[1](i, j+h, k) = this->dt / this->dj * 
                                        (
                                         10 * a / tau * sin(2 * xp) * cos(y1) * cos(pi * t / tau)
                                        ) * cos(y1) * rho;
          
          this->mem->GC[2](i, j, k+h) = this->dt / this->dk * 
                                        (
                                         omg_0 * sin(xp) * cos(y) * cos(2 * pi * t / tau) * s1
                                        ) * a * cos(y) * (-1.0) / g;


          if (needs_dt)
          {
            this->mem->ndt_GC[0](i+h, j, k) = this->dt / this->di * 
                                        (
                                         - 10 * a / tau * pow(sin(xp1), 2) * sin(2 * y) * pi / tau * sin(pi * t / tau)
                                         - udv * 2 * pi / tau * sin(2 * pi * t / tau)
                                        ) * rho * this->dt;

            this->mem->ndt_GC[1](i, j+h, k) = this->dt / this->dj *  
                                        (
                                         -10 * a / tau * sin(2 * xp) * cos(y1) * pi / tau * sin(pi * t / tau)
                                        ) * cos(y1) * rho * this->dt;

            this->mem->ndt_GC[2](i, j, k+h) = this->dt / this->dk *  
                                        (
                                         -omg_0 * sin(xp) * cos(y) * 2 * pi / tau * sin(2 * pi * t / tau) * s1
                                        ) * a * cos(y) * (-1.0) / g * this->dt;
            
            this->mem->ndtt_GC[0](i+h, j, k) = this->dt / this->di * 
                                        (
                                         - 10 * a / tau * pow(sin(xp1), 2) * sin(2 * y) * pow(pi / tau, 2) * cos(pi * t / tau)
                                         - udv * pow(2 * pi / tau, 2) * cos(2 * pi * t / tau)
                                        ) * rho * pow(this->dt, 2);

            this->mem->ndtt_GC[1](i, j+h, k) = this->dt / this->dj *  
                                        (
                                         -10 * a / tau * sin(2 * xp) * cos(y1) * pow(pi / tau, 2) * cos(pi * t / tau)
                                        ) * cos(y1) * rho * pow(this->dt, 2);

            this->mem->ndtt_GC[2](i, j, k+h) = this->dt / this->dk *  
                                        (
                                         -omg_0 * sin(xp) * cos(y) * pow(2 * pi / tau, 2) * cos(2 * pi * t / tau) * s1
                                        ) * a * cos(y) * (-1.0) / g * pow(this->dt, 2);
          }
        }
      }
    }

    auto ex = this->halo - 1;
    this->xchng_vctr_alng(this->mem->GC);
    this->xchng_vctr_nrml(this->mem->GC, this->i^ex, this->j^ex, this->k^ex);

    if (needs_dt)
    {
      this->xchng_vctr_alng(this->mem->ndt_GC);
      this->xchng_vctr_nrml(this->mem->ndt_GC, this->i^ex, this->j^ex, this->k^ex);
      
      this->xchng_vctr_alng(this->mem->ndtt_GC);
      this->xchng_vctr_nrml(this->mem->ndtt_GC, this->i^ex, this->j^ex, this->k^ex);
    }

    return true;
  }

  public:
  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t pi, a, g, z_top, p_top, p_0, rho_0, tau, H, omg_0, b;
  };

  // ctor
  dcmip_1( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p),
    pi(p.pi),
    a(p.a),
    g(p.g),
    z_top(p.z_top),
    p_top(p.p_top),
    p_0(p.p_0),
    rho_0(p.rho_0),
    tau(p.tau),
    H(p.H),
    omg_0(p.omg_0),
    b(p.b)
  {}
};
