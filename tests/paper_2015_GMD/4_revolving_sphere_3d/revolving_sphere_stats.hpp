/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * errors for over the sphere test
 */
#pragma once

template<class parent_t>
struct stats : public parent_t
{
  using parent_t::parent_t;
  using real_t = typename parent_t::real_t;

  typename parent_t::arr_t ic, fc;

  real_t last_timestep;

  std::ofstream ofs;

  void hook_ante_loop(const int nt)
  {
    parent_t::hook_ante_loop(nt);
    if (this->rank != 0) return;

    //checking what are the MPDATA options of each test simulation (iga / fct / ...) 
    //basing on hdf outdir name and naming output stats file accordingly
    if(!ofs.is_open())
      ofs.open("stats_" + this->outdir + ".txt", std::ofstream::out);

    last_timestep = nt;

    ic.resize(this->mem->advectee_global().shape());
    ic = this->mem->advectee_global();

    // at the end of the test true solution is equal to the inital state
    ofs << std::fixed << std::setprecision(8) << std::endl;
    ofs << "timestep      = 0"                << std::endl;
    ofs << "min(solution) = " << min(ic)      << std::endl;
    ofs << "max(solution) = " << max(ic)      << std::endl;
    ofs << " "                                << std::endl;
  }

  void hook_post_step()
  {
    parent_t::hook_post_step();
    this->mem->barrier();
    if (this->rank != 0) return;
    if (this->timestep == last_timestep) 
    { 
      // final condition
      fc.resize(this->mem->advectee_global().shape());
      fc = this->mem->advectee_global();

      //output stats
      ofs << "timestep = " << this->timestep                                                                                  << std::endl;
      ofs << "max(solution) = " << max(fc)                                                                                    << std::endl;
      ofs << "min(solution) = " << min(fc)                                                                                    << std::endl;
      ofs << "Linf = " << max(abs(fc - ic))                                           /*nx */                                 << std::endl;
      ofs << "L2 = " << 1.0 / (this->timestep * this->dt) * sqrt(1.0 / pow(this->mem->distmem.grid_size[0], 3) * sum(pow2(fc - ic)))  << std::endl;
      ofs << "L1 = " << 1.0 / (this->timestep * this->dt) *     (1.0 / pow(this->mem->distmem.grid_size[0], 3) * sum(abs (fc - ic)))  << std::endl; 
    }
  }
};

