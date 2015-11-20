/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <cmath>

template <class ct_params_t>
class pbl : public libmpdataxx::solvers::mpdata_rhs_vip_prs<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata_rhs_vip_prs<ct_params_t>;
  using ix = typename ct_params_t::ix;

  public:
  using real_t = typename ct_params_t::real_t;

  private:
  typename parent_t::arr_t &Tht_e, &H, &tmp1, &tmp2;
  real_t g, Tht_ref;

  void update_rhs(
    libmpdataxx::arrvec_t<
      typename parent_t::arr_t
    > &rhs, 
    const real_t &dt, 
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at); 

    const auto &Tht = this->state(ix::tht); 
    const auto &u = this->state(ix::u); 
    const auto &v = this->state(ix::v); 
    const auto &w = this->state(ix::w); 
    const auto &ijk = this->ijk;
    const auto &i = this->i;
    const auto &j = this->j;
    const auto &k = this->k;

    switch (at) 
    {
      case (0): 
      {
        rhs.at(ix::tht)(ijk) += H(ijk) - (*this->mem->vab_coeff)(ijk) * (Tht(ijk) - Tht_e(ijk));

        tmp1(ijk) = g * (Tht(ijk) - Tht_e(ijk)) / Tht_ref;
        this->xchng_sclr(tmp1, i, j, k);
        tmp2(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j ,k) + tmp1(i, j, k - 1));
        rhs.at(ix::w)(ijk) += tmp2(ijk);
        break;
      }
   
      case (1): 
      {
        rhs.at(ix::tht)(ijk) += (H(ijk) - (*this->mem->vab_coeff)(ijk) * (Tht(ijk) - Tht_e(ijk)))
                                / (1.0 + 0.5 * (*this->mem->vab_coeff)(ijk) * this->dt);

        tmp1(ijk) = g * (Tht(ijk) + 0.5 * this->dt * rhs.at(ix::tht)(ijk) - Tht_e(ijk)) / Tht_ref;
        this->xchng_sclr(tmp1, i, j, k);
        tmp2(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j ,k) + tmp1(i, j, k - 1));
        rhs.at(ix::w)(ijk) += tmp2(ijk);
        break;
      }
    }

  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t g = 10.0, Tht_ref = 0; 
    typename parent_t::arr_t *Tht_e, *H;
  };

  // ctor
  pbl( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p),
    Tht_e(*p.Tht_e),
    H(*p.H),
    g(p.g),
    Tht_ref(p.Tht_ref),
    tmp1(args.mem->tmp[__FILE__][0][0]),
    tmp2(args.mem->tmp[__FILE__][0][1])
  {
    assert(Tht_ref != 0);
  }

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // tmp1, tmp2
  }
};
