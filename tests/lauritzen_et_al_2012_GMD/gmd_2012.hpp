/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include <libmpdata++/solvers/mpdata.hpp>

using T = double;
constexpr T pi = boost::math::constants::pi<T>();

template <class ct_params_t, class ct_test_params_t>
class gmd_2012 : public libmpdataxx::solvers::mpdata<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata<ct_params_t>;
  using parent_t::parent_t;
  
  public:

  using real_t = typename ct_params_t::real_t;

  protected:

  using tp = ct_test_params_t;

  libmpdataxx::arrvec_t<typename parent_t::arr_t> &gc_node;

  std::pair<real_t, real_t> gc_exact(const real_t t, const real_t x, const real_t y)
  {
    const auto xp = x - 2 * pi * t / tp::T;

    return {
             (   10 * tp::R / tp::T * pow(sin(xp), 2) * sin(2 * y) * cos(pi * t / tp::T)
               + 2 * pi * tp::R / tp::T * cos(y)
             ) * this->dj * this->dt
           ,
             (   10 * tp::R / tp::T * sin(2 * xp) * pow(cos(y), 2) * cos(pi * t / tp::T)
             ) * this->di * this->dt
           };
  }

  bool calc_gc() final
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

        auto gc_ij = gc_exact(t, x, y);

        gc_node[0](i, j) = gc_ij.first;
        gc_node[1](i, j) = gc_ij.second;
      }
    }

    //this->mem->barrier();
    this->xchng_sclr(gc_node[0], this->i, this->j);
    this->xchng_sclr(gc_node[1], this->i, this->j);

    // calculate staggered advectors
    for (int i = this->i.first()-1; i <= this->i.last(); ++i)
    {
      for (int j = this->j.first()-1; j <= this->j.last(); ++j)
      {
        //const auto x = i * this->di;
        //const auto y = (j+ 0.5) * this->dj - pi / 2;
        //auto gc_ij_x = gc_exact(t, x + 0.5 * this->di, y);
        //auto gc_ij_y = gc_exact(t, x, y + 0.5 * this->dj);
        //this->mem->GC[0](i+h, j) = gc_ij_x.first;
        //this->mem->GC[1](i, j+h) = gc_ij_y.second;

        this->mem->GC[0](i+h, j) = 0.5 * (gc_node[0](i, j) + gc_node[0](i+1, j));
        this->mem->GC[1](i, j+h) = 0.5 * (gc_node[1](i, j) + gc_node[1](i, j+1));
      }
    }

    auto ex = this->halo - 1;
    //this->xchng_vctr_alng(this->mem->GC);
    this->xchng_vctr_nrml(this->mem->GC, this->i^ex, this->j^ex);
    
    //auto div = max(abs(
    //                   ( this->mem->GC[0](this->i+h, this->j) - this->mem->GC[0](this->i-h, this->j)  
    //                   + this->mem->GC[1](this->i, this->j+h) - this->mem->GC[1](this->i, this->j-h)
    //                   ) / (*this->mem->G)(this->i, this->j)
    //                  ));

    //if (this->time == 0)
    //{
    //  std::stringstream ss;
    //  ss << "rank: " << this->rank << " div: " << div << std::endl;
    //  std::cout << ss.str() << std::endl;
    //}

    return true;
  }

  public:

  // ctor
  gmd_2012( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p),
    gc_node(args.mem->tmp[__FILE__][0])
  {}

  static void alloc(
    typename parent_t::mem_t *mem, 
    const int &n_iters
  ) {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // gc_node
  }
};
