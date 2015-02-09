/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * output of errors and extrema for boussinesq test
 */
#pragma once
#include <fstream>

template <class parent_t>
struct stats : public parent_t
{
  using parent_t::parent_t;
  using real_t = typename parent_t::real_t;
  using ix = typename parent_t::ix;

  real_t 
    sum_init_pert,
    sum_init_pert_2,
    sum_init_tht,
    sum_init_tht_2;

  std::ofstream ofs;

  void hook_ante_loop(const int nt)
  {
    parent_t::hook_ante_loop(nt);
    if (this->mem->rank() != 0) return;

    if (!ofs.is_open())
      ofs.open("stats.txt", std::ofstream::out);

    sum_init_pert   = sum(// due to cyclic boundary condition summing has to be done over the whole domain minus one edge per each dimension
      (this->mem->advectee(ix::tht) - this->Tht_ref)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd))
    );
    sum_init_pert_2 = sum(// due to cyclic boundary condition summing has to be done over the whole domain minus one edge per each dimension
      pow((this->mem->advectee(ix::tht) - this->Tht_ref)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd)) , 2)
    );
    sum_init_tht    = sum(// due to cyclic boundary condition summing has to be done over the whole domain minus one edge per each dimension
      this->mem->advectee(ix::tht)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd))
    );
    sum_init_tht_2  = sum(// due to cyclic boundary condition summing has to be done over the whole domain minus one edge per each dimension
      pow(this->mem->advectee(ix::tht)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd)), 2)
    );

    // outputting definition of errors ...
    ofs << "error1 := (sum of perturbation - sum of initial perturbation) / sum of initial perturbation * 100%" << std::endl;
    ofs << "error2 := (sum of perturbation^2 - sum of initial perturbation^2) / sum of initial perturbation^2 * 100%" << std::endl;
    ofs << "error3 := (sum of theta - sum of initial theta) / sum of initial theta * 100%" << std::endl;
    ofs << "error4 := (sum of theta^2 - sum of initial theta^2) / sum of initial theta^2 * 100%" << std::endl;
    ofs << " " << std::endl;     
    // ... initial perturbation ...
    ofs << std::scientific << std::setprecision(6)<<std::endl;
    ofs <<"timestep = 0" << std::endl;
    ofs << "sum of initial perturbation         = " << sum_init_pert   << std::endl;
    ofs << "sum of initial perturbation squared = " << sum_init_pert_2 << std::endl;
    ofs << "sum of initial theta                = " << sum_init_tht    << std::endl;
    ofs << "sum of initial theta squared        = " << sum_init_tht_2  << std::endl;
    ofs << " " << std::endl; 
  }

  void hook_post_step()
  {
    parent_t::hook_post_step();
    this->mem->barrier();
    if (this->mem->rank() != 0) return;

    if (this->timestep % 100 == 0)
    {
      if (!ofs.is_open())
        ofs.open("stats.txt", std::ofstream::out);
      // ... and errors and signal extrema throughout the simulation
//      ofs << std::fixed << std::setprecision(8)<<std::endl;
      ofs << "timestep      = " << this->timestep<<std::endl;
      ofs << "min(tht)      = " << min(this->mem->advectee(ix::tht)) << std::endl;
      ofs << "max(tht)      = " << max(this->mem->advectee(ix::tht)) << std::endl;
      ofs << "max(velocity) = " << max(sqrt(pow(this->mem->advectee(ix::u),2) + pow(this->mem->advectee(ix::w),2))) <<std::endl;
      ofs << "error1        = " << 
        (                                                    //see ante_loop comment 
          sum((this->mem->advectee(ix::tht) - this->Tht_ref)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd)))  - sum_init_pert
        ) / sum_init_pert * 100 << std::endl;
      ofs << "error2        = " << 
        (
          sum(
            pow((this->mem->advectee(ix::tht) - this->Tht_ref)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd)), 2)
          )  - sum_init_pert_2
        ) / sum_init_pert_2 * 100 << std::endl;
      ofs << "error3        = " << 
        (
          sum(this->mem->advectee(ix::tht)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd)))  - sum_init_tht
        ) / sum_init_tht * 100 << std::endl;
      ofs << "error4        = " << 
        (
          sum(
            pow(this->mem->advectee(ix::tht)(blitz::Range(1, blitz::toEnd), blitz::Range(1, blitz::toEnd)),2)
          )  - sum_init_tht_2
        ) / sum_init_tht_2 * 100 << std::endl;
      ofs << "  " << std::endl; 
    }
  }
};
