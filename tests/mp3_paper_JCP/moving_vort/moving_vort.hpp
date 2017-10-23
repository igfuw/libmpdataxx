#pragma once
#include <libmpdata++/solvers/mpdata.hpp>
#include "transforms.hpp"

template <class ct_params_t, class ct_test_params_t>
class moving_vort : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  using tp = ct_test_params_t;

  static constexpr bool needs_dt = libmpdataxx::opts::isset(ct_params_t::opts, libmpdataxx::opts::div_3rd) ||
                                   libmpdataxx::opts::isset(ct_params_t::opts, libmpdataxx::opts::div_3rd_dt);

  real_t xc, yc;

  libmpdataxx::arrvec_t<typename parent_t::arr_t> &ex_node;

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

    if (needs_dt)
    {
      update_vortex_centre(t - 0.5 * this->dt);
      for (int i = this->i.first(); i <= this->i.last(); ++i)
      {
        for (int j = this->j.first(); j <= this->j.last(); ++j)
        {
          const auto x = i * this->di;
          const auto y = (j+ 0.5) * this->dj - pi / 2;

          auto gc_ij = gc_exact(x, y);

          ex_node[2](i, j) = -gc_ij.first;
          ex_node[3](i, j) = -gc_ij.second;
          ex_node[4](i, j) = gc_ij.first;
          ex_node[5](i, j) = gc_ij.second;
        }
      }
      
      update_vortex_centre(t);
      for (int i = this->i.first(); i <= this->i.last(); ++i)
      {
        for (int j = this->j.first(); j <= this->j.last(); ++j)
        {
          const auto x = i * this->di;
          const auto y = (j+ 0.5) * this->dj - pi / 2;

          auto gc_ij = gc_exact(x, y);

          ex_node[4](i, j) -= 2 * gc_ij.first;
          ex_node[5](i, j) -= 2 * gc_ij.second;
        }
      }
    }

    update_vortex_centre(t + 0.5 * this->dt);

    // calculate exact advectors at grid points
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        const auto x = i * this->di;
        const auto y = (j+ 0.5) * this->dj - pi / 2;

        auto gc_ij = gc_exact(x, y);

        ex_node[0](i, j) = gc_ij.first;
        ex_node[1](i, j) = gc_ij.second;

        if (needs_dt)
        {
          ex_node[2](i, j) += gc_ij.first;
          ex_node[3](i, j) += gc_ij.second;

          ex_node[4](i, j) += gc_ij.first;
          ex_node[5](i, j) += gc_ij.second;
          ex_node[4](i, j) *= 4;
          ex_node[5](i, j) *= 4;
        }
      }
    }

    for (int d = 0; d < 2; ++d)
    {
      this->xchng_sclr(ex_node[d], this->ijk);
      if (needs_dt)
      {
        this->xchng_sclr(ex_node[2 + d], this->ijk);
        this->xchng_sclr(ex_node[4 + d], this->ijk);
      }
    }

    // calculate staggered advectors
    auto ex = this->halo - 1;
    for (int i = this->i.first()-1-ex; i <= this->i.last()+ex; ++i)
    {
      for (int j = this->j.first()-1-ex; j <= this->j.last()+ex; ++j)
      {
        this->mem->GC[0](i+h, j) = 0.5 * (ex_node[0](i, j) + ex_node[0](i+1, j));
        this->mem->GC[1](i, j+h) = 0.5 * (ex_node[1](i, j) + ex_node[1](i, j+1));
        if (needs_dt)
        {
          this->mem->ndt_GC[0](i+h, j) = 0.5 * (ex_node[2](i, j) + ex_node[2](i+1, j));
          this->mem->ndt_GC[1](i, j+h) = 0.5 * (ex_node[3](i, j) + ex_node[3](i, j+1));

          this->mem->ndtt_GC[0](i+h, j) = 0.5 * (ex_node[4](i, j) + ex_node[4](i+1, j));
          this->mem->ndtt_GC[1](i, j+h) = 0.5 * (ex_node[5](i, j) + ex_node[5](i, j+1));
        }
      }
    }

    this->xchng_vctr_nrml(this->mem->GC, this->ijk, ex);
    if (needs_dt)
    {
      this->xchng_vctr_nrml(this->mem->ndt_GC, this->ijk, ex);
      this->xchng_vctr_nrml(this->mem->ndtt_GC, this->ijk, ex);
    }

    return true;
  }

  public:

  // ctor
  moving_vort( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p),
    ex_node(args.mem->tmp[__FILE__][0]),
    xpf{tp::x0, tp::y0},
    ypf{tp::x0, tp::y0}
  {}

  static void alloc(
    typename parent_t::mem_t *mem, 
    const int &n_iters
  ) {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, needs_dt ? 6 : 2); // ex_node
  }
};
