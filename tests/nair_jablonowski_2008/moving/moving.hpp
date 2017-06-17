/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>
#include "../common/transforms.hpp"

template <class ct_params_t, class ct_test_params_t>
class moving : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  using tp = ct_test_params_t;

  real_t xc, yc;

  libmpdataxx::arrvec_t<typename parent_t::arr_t> &gc_node;

  xpf_t xpf;
  ixpf_t ixpf;
  ypf_t ypf;
  iypf_t iypf;

  void update_vortex_centre(real_t t)
  {
    xpf.x0 = pi;
    xpf.y0 = pi / 2 - tp::a;
    ypf.x0 = pi;
    ypf.y0 = pi / 2 - tp::a;

    ixpf.x0 = pi;
    ixpf.y0 = pi / 2 - tp::a;
    iypf.x0 = pi;
    iypf.y0 = pi / 2 - tp::a;

    real_t x = xpf(tp::x0, tp::y0);
    real_t y = ypf(tp::x0, tp::y0);

    x += tp::u0 * t;

    xc = ixpf(x, y);
    yc = iypf(x, y);
    
    xpf.x0 = xc;
    xpf.y0 = yc;
    ypf.x0 = xc;
    ypf.y0 = yc;
  }

  std::pair<real_t, real_t> gc_exact(const real_t x, const real_t y)
  {
    const real_t r = 3 * cos(ypf(x, y));
    const real_t omg = r != 0 ? tp::v0 * 3 * sqrt(2.) / (2 * r) * tanh(r) / pow(cosh(r), 2) : 0;
    return {
             (   omg * (sin(yc) - cos(yc) * cos(x - xc) * tan(y))
             + tp::u0 * (cos(tp::a) + tan(y) * sin(tp::a) * cos(x))
             ) * cos(y) * this->dj * this->dt
           ,
             (   omg * cos(yc) * sin(x - xc)
             - tp::u0 * sin(x) * sin(tp::a)
             ) * cos(y) * this->di * this->dt
           };
  }

  bool calc_gc() final
  {
    using namespace libmpdataxx::arakawa_c;

    auto t = this->time;
    update_vortex_centre(t + 0.5 * this->dt);

    // calculate exact advectors at grid points
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        const auto x = i * this->di;
        const auto y = (j+ 0.5) * this->dj - pi / 2;

        auto gc_ij = gc_exact(x, y);

        gc_node[0](i, j) = gc_ij.first;
        gc_node[1](i, j) = gc_ij.second;
      }
    }

    this->mem->barrier();

    // calculate staggered advectors
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        this->mem->GC[0](i+h, j) = 0.5 * (gc_node[0](i, j) + gc_node[0](i+1, j));
        this->mem->GC[1](i, j+h) = 0.5 * (gc_node[1](i, j) + gc_node[1](i, j+1));
      }
    }

    auto ex = this->halo - 1;
    this->xchng_vctr_alng(this->mem->GC);
    this->xchng_vctr_nrml(this->mem->GC, this->ijk, ex);

    return true;
  }

  public:

  // ctor
  moving( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p),
    gc_node(args.mem->tmp[__FILE__][0]),
    xpf{tp::x0, tp::y0},
    ypf{tp::x0, tp::y0}
  {}

  static void alloc(
    typename parent_t::mem_t *mem, 
    const int &n_iters
  ) {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // gc_node
  }
};
