#pragma once
#include <libmpdata++/solvers/mpdata.hpp>

using T = double;
constexpr T pi = boost::math::constants::pi<T>();

template <class ct_params_t, class ct_test_params_t>
class reversing_deform : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  using tp = ct_test_params_t;

  static constexpr bool needs_dt = libmpdataxx::opts::isset(ct_params_t::opts, libmpdataxx::opts::div_3rd) ||
                                   libmpdataxx::opts::isset(ct_params_t::opts, libmpdataxx::opts::div_3rd_dt);

  libmpdataxx::arrvec_t<typename parent_t::arr_t> &ex_node;

  struct exact_t
  {
    real_t gc[2];
    real_t ndt_gc[2];
    real_t ndtt_gc[2];
  };

  exact_t exact(const real_t t, const real_t x1, const real_t y1, const real_t x2, const real_t y2)
  {
    using std::sin;
    using std::cos;
    using std::pow;
    exact_t res;

    const auto x1p = x1 - 2 * pi * t / tp::T;
    const auto x2p = x2 - 2 * pi * t / tp::T;

    res.gc[0] = (   10 * tp::R / tp::T * pow(sin(x1p), 2) * sin(2 * y1) * cos(pi * t / tp::T)
                  + 2 * pi * tp::R / tp::T * cos(y1)
                ) * this->dj * this->dt;
    res.gc[1] = (   10 * tp::R / tp::T * sin(2 * x2p) * pow(cos(y2), 2) * cos(pi * t / tp::T)
                ) * this->di * this->dt;

    if (needs_dt)
    {
      res.ndt_gc[0] = (   -10 * tp::R * pi / pow(tp::T, 2) * pow(sin(x1p), 2) * sin(2 * y1) * sin(pi * t / tp::T)
                      ) * this->dj * pow(this->dt, 2);
      res.ndt_gc[1] = (   -10 * tp::R * pi / pow(tp::T, 2) * sin(2 * x2p) * pow(cos(y2), 2) * sin(pi * t / tp::T)
                      ) * this->di * pow(this->dt, 2);
      
      res.ndtt_gc[0] = (   -10 * tp::R * pow(pi, 2) / pow(tp::T, 3) * pow(sin(x1p), 2) * sin(2 * y1) * cos(pi * t / tp::T)
                       ) * this->dj * pow(this->dt, 3);
      res.ndtt_gc[1] = (   -10 * tp::R * pow(pi, 2) / pow(tp::T, 3) * sin(2 * x2p) * pow(cos(y2), 2) * cos(pi * t / tp::T)
                       ) * this->di * pow(this->dt, 3);
    }
    return res;
  }

  void calc_gc_staggered()
  {
    using namespace libmpdataxx::arakawa_c;
    const auto t = this->time + 0.5 * this->dt;

    // calculate staggered advectors
    for (int i = this->i.first()-1; i <= this->i.last(); ++i)
    {
      for (int j = this->j.first()-1; j <= this->j.last(); ++j)
      {
        const auto x1 = (i + 0.5) * this->di;
        const auto y1 = (j+ 0.5) * this->dj - pi / 2;
        
        const auto x2 = i * this->di;
        const auto y2 = (j+ 1.0) * this->dj - pi / 2;

        const auto ex = exact(t, x1, y1, x2, y2);
        
        this->mem->GC[0](i+h, j) = ex.gc[0];
        this->mem->GC[1](i, j+h) = ex.gc[1];
        if (needs_dt)
        {
          this->mem->ndt_GC[0](i+h, j) = ex.ndt_gc[0];
          this->mem->ndt_GC[1](i, j+h) = ex.ndt_gc[1];

          this->mem->ndtt_GC[0](i+h, j) = ex.ndtt_gc[0];
          this->mem->ndtt_GC[1](i, j+h) = ex.ndtt_gc[1];
        }
      }
    }
    
    auto ex = this->halo - 1;
    this->xchng_vctr_alng(this->mem->GC);
    this->xchng_vctr_nrml(this->mem->GC, this->ijk, ex);
    if (needs_dt)
    {
      this->xchng_vctr_alng(this->mem->ndt_GC);
      this->xchng_vctr_alng(this->mem->ndtt_GC);
      this->xchng_vctr_nrml(this->mem->ndt_GC, this->ijk, ex);
      this->xchng_vctr_nrml(this->mem->ndtt_GC, this->ijk, ex);
    }
  }
  
  void calc_gc_node()
  {
    using namespace libmpdataxx::arakawa_c;
    const auto t = this->time + 0.5 * this->dt;

    // calculate exact advectors at grid points
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        const auto x = i * this->di;
        const auto y = (j+ 0.5) * this->dj - pi / 2;

        const auto ex = exact(t, x, y, x, y);
        
        for (int d = 0; d < 2; ++d)
        {
          ex_node[d](i, j) = ex.gc[d];
        
          if (needs_dt)
          {
            ex_node[2 + d](i, j) = ex.ndt_gc[d];
            ex_node[4 + d](i, j) = ex.ndtt_gc[d];
          }
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
  }

  bool calc_gc() final
  {
    if (tp::node) 
      calc_gc_node();
    else 
      calc_gc_staggered();
    return true;
  }

  public:

  // ctor
  reversing_deform( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p),
    ex_node(args.mem->tmp[__FILE__][0])
  {}

  static void alloc(
    typename parent_t::mem_t *mem, 
    const int &n_iters
  ) {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, needs_dt ? 6 : 2); // ex_node
  }
};
