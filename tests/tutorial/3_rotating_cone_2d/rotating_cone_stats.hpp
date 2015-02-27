/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * errors and signal extrema for rotating cone 2D test
 */
#pragma once

template<class parent_t>
struct stats : public parent_t
{
  using parent_t::parent_t;
  using real_t = typename parent_t::real_t;

  typename parent_t::arr_t true_solution, error_2;

  real_t last_timestep;

  std::ofstream ofs;
 
  void hook_ante_loop(const int nt)
  {
    parent_t::hook_ante_loop(nt);
    if (this->mem->rank() != 0) return;

    //checking what are the MPDATA options of each test simulation (fct / iga / ...) 
    //basing on gnuplot output filename ...
    std::string gnuplot_name = this->p.gnuplot_output;
    std::string name;
    char delimeter('%');
    std::istringstream iss(gnuplot_name);
    getline(iss, name, delimeter);
    // ... and naming output stats file accordingly
    if(!ofs.is_open())
      ofs.open("stats_"+name+".txt", std::ofstream::out);

    last_timestep = nt;

    true_solution.resize(this->mem->advectee().shape());
    error_2.resize(this->mem->advectee().shape());

    //after one full rotation true solution is equal to the inital state
    true_solution = this->mem->advectee();
    ofs << std::fixed << std::setprecision(8) << std::endl;
    ofs << "timestep      = 0"     << std::endl;
    ofs << "min(solution) = " << min(true_solution) <<std::endl;
    ofs << "max(solution) = " << max(true_solution) <<std::endl;
    ofs << " " <<std::endl;
  }

  void hook_post_step()
  {
    parent_t::hook_post_step();
    this->mem->barrier();
    if (this->mem->rank() != 0) return;
    if (this->timestep == last_timestep) 
    { 
      ofs << "timestep     = " << this->timestep << std::endl;
      ofs << "min(psi)     = " << min(this->mem->advectee()) << std::endl;
      ofs << "max(psi)     = " << max(this->mem->advectee()) << std::endl;
 
      error_2 = pow(true_solution - this->mem->advectee(), 2);
      ofs << "rms error    = " 
        << sqrt(sum(error_2) / this->mem->advectee().extent(0) / this->mem->advectee().extent(1)) / this->timestep / this->dt << std::endl;
    }
  }
};
