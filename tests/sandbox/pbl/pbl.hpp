/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/boussinesq.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include <libmpdata++/formulae/refined_grid.hpp>
#include <cmath>

template <class ct_params_t>
class pbl : public libmpdataxx::output::hdf5_xdmf<libmpdataxx::solvers::boussinesq<ct_params_t>>
{
  using ix = typename ct_params_t::ix;

  public:
  using real_t = typename ct_params_t::real_t;
  using parent_t = libmpdataxx::output::hdf5_xdmf<libmpdataxx::solvers::boussinesq<ct_params_t>>;

  private:
  real_t hscale, iles_cdrag;
  typename parent_t::arr_t &tke, &r2r_avg;

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

      // output tht refined with fractal reconstruction
      this->generate_stretching_parameters(std::random_device{}(), libmpdataxx::formulae::fractal::stretch_params::d_distro_t::DNS_vel);
      this->reconstruct_refinee(ix::w);

//      this->generate_stretching_parameters(std::random_device{}(), libmpdataxx::formulae::fractal::stretch_params::d_distro_t::LES_rv_supersaturated);
//      this->generate_stretching_parameters(std::random_device{}(), libmpdataxx::formulae::fractal::stretch_params::d_distro_t::LES_th_supersaturated);
//      this->reconstruct_refinee(ix::tht);

      this->mem->barrier();
      if (this->rank == 0)
      {
        if (static_cast<libmpdataxx::solvers::sgs_scheme_t>(ct_params_t::sgs_scheme) == libmpdataxx::solvers::smg)
        {
          this->record_aux_dsc("tke", this->tke);
        }
        this->record_aux_dsc("p", this->Phi);
        this->record_aux_dsc_refined("w reconstructed using DNS_vel", this->mem->refinee(this->ix_r2r.at(ix::w)));
//        this->record_aux_dsc_refined("tht reconstructed", this->mem->refinee(this->ix_r2r.at(ix::tht)));
      }
      this->mem->barrier();

      // another reconstruction of w using different stretchin parameteres

      this->generate_stretching_parameters(std::random_device{}(), libmpdataxx::formulae::fractal::stretch_params::d_distro_t::LES_th_supersaturated);
      this->reconstruct_refinee(ix::w);

      this->mem->barrier();
      if (this->rank == 0)
      {
        this->record_aux_dsc_refined("w reconstructed using LES_th_supersaturated", this->mem->refinee(this->ix_r2r.at(ix::w)));
      }
      this->mem->barrier();

      // output tht refined with linear interpolation
      this->interpolate_refinee(ix::tht);

      this->mem->barrier();
      if (this->rank == 0)
      {
        this->record_aux_dsc_refined("tht interpolated", this->mem->refinee(this->ix_r2r.at(ix::tht)));
      }
      this->mem->barrier();

      // output tht on refined grid after averaging from interpolated refined tht to regular grid tht
      // TODO: fill refined grid before this average!
      // TODO: this will make avg_edge_sclr obsolete?
      libmpdataxx::formulae::refined::spatial_average_ref2reg<real_t>(this->mem->refinee(this->ix_r2r.at(ix::tht)), this->ijk_r2r, this->mem->n_ref/2, this->mem->distmem.grid_size_ref, true);
      this->r2r_avg(this->ijk) = this->mem->refinee(this->ix_r2r.at(ix::tht))(this->ijk_r2r); 

      this->mem->barrier();
      if (this->rank == 0)
      {
        //this->record_aux_dsc_refined("tht interpolated and averaged", this->mem->refinee(this->ix_r2r.at(ix::tht)));
        this->record_aux_dsc("tht averaged from interpolated refined grid", this->r2r_avg);
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
    tke(args.mem->tmp[__FILE__][0][0]),
    r2r_avg(args.mem->tmp[__FILE__][0][1])
  {}

  static void alloc(
    typename parent_t::mem_t *mem, 
    const int &n_iters
  ) {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // tke
  }
};
