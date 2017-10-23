/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/boussinesq.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include <cmath>

template <class ct_params_t>
class pbl : public libmpdataxx::output::hdf5_xdmf<libmpdataxx::solvers::boussinesq<ct_params_t>>
{
  using parent_t = libmpdataxx::output::hdf5_xdmf<libmpdataxx::solvers::boussinesq<ct_params_t>>;
  using ix = typename ct_params_t::ix;

  public:
  using real_t = typename ct_params_t::real_t;

  private:
  real_t hscale, iles_cdrag;
  typename parent_t::arr_t &tke;

  void multiply_sgs_visc()
  {
    parent_t::multiply_sgs_visc();

    if (this->timestep % static_cast<int>(this->outfreq) == 0 &&
        static_cast<libmpdataxx::solvers::sgs_scheme_t>(ct_params_t::sgs_scheme) == libmpdataxx::solvers::smg)
    {
      tke(this->ijk) = pow2(this->k_m(this->ijk) / (this->c_m * this->mix_len(this->ijk)));
    }
  }

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();

    if (static_cast<libmpdataxx::solvers::sgs_scheme_t>(ct_params_t::sgs_scheme) == libmpdataxx::solvers::iles)
    {
      for (int k = this->k.first(); k <= this->k.last(); ++k)
      {
        this->vip_rhs[0](this->i, this->j, k) += - 2 * iles_cdrag / hscale * sqrt(
                                                  pow2(this->state(ix::vip_i)(this->i, this->j, 0))
                                                + pow2(this->state(ix::vip_j)(this->i, this->j, 0))
                                                ) * this->state(ix::vip_i)(this->i, this->j, 0)
                                                  * exp(-this->dj * k / hscale);
        
        this->vip_rhs[1](this->i, this->j, k) += - 2 * iles_cdrag / hscale * sqrt(
                                                  pow2(this->state(ix::vip_i)(this->i, this->j, 0))
                                                + pow2(this->state(ix::vip_j)(this->i, this->j, 0))
                                                ) * this->state(ix::vip_j)(this->i, this->j, 0)
                                                  * exp(-this->dj * k / hscale);
      }
    }
    
    if (this->timestep % static_cast<int>(this->outfreq) == 0)
    {
      if (this->rank == 0) std::cout << this->timestep << std::endl;
      this->mem->barrier();
      if (this->rank == 0)
      {
        if (static_cast<libmpdataxx::solvers::sgs_scheme_t>(ct_params_t::sgs_scheme) == libmpdataxx::solvers::smg)
        {
          this->record_aux_dsc("tke", this->tke);
        }
        this->record_aux_dsc("p", this->Phi);
      }
      this->mem->barrier();
    }
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t hscale = 1, iles_cdrag = 0; 
  };

  // ctor
  pbl( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p),
    hscale(p.hscale),
    iles_cdrag(p.iles_cdrag),
    tke(args.mem->tmp[__FILE__][0][0])
  {}

  static void alloc(
    typename parent_t::mem_t *mem, 
    const int &n_iters
  ) {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 1); // tke
  }
};
