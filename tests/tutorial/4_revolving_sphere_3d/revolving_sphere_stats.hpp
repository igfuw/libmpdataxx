/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * errors for over the sphere test
 */
#pragma once
  // initial condition
  auto ic = slv.advectee().copy();

  // time stepping
  slv.advance(nt);

 // final condition
  auto fc = slv.advectee();

  std::cout << "max:\t" << max(fc) << std::endl;
  std::cout << "min:\t" << min(fc) << std::endl;
  std::cout << "Linf:\t" << max(abs(fc - ic)) << std::endl;
  std::cout << "L2:\t" << 1.0 / (nt * dt) * sqrt(1.0 / pow(nx, 3) * sum(pow2(fc - ic))) << std::endl;
  std::cout << "L1:\t" << 1.0 / (nt * dt) * (1.0 / pow(nx, 3) * sum(abs(fc - ic))) << std::endl; 

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
    if (this->mem->rank() != 0) return;

    //checking what are the MPDATA options of each test simulation (default / best) 
    //basing on hdf outdir name and naming output stats file accordingly
    if(!ofs.is_open())
      ofs.open("stats_" + this->outdir + ".txt", std::ofstream::out);

    last_timestep = nt;

    ic.resize(ii, jj);
    ic = this->mem->advectee().copy()(ii, jj);

    // at the end of the test true solution is equal to the inital state
    ofs << std::fixed << std::setprecision(8) << std::endl;
    ofs << "timestep      = 0"     << std::endl;
    ofs << "min(solution) = " << min(ic) <<std::endl;
    ofs << "max(solution) = " << max(ic) <<std::endl;
    ofs << " " <<std::endl;
  }

  void hook_post_step()
  {
    parent_t::hook_post_step();
    this->mem->barrier();
    if (this->mem->rank() != 0) return;
    if (this->timestep == last_timestep) 
    { 
      // final condition
      fc.resize(ii, jj);
      fc = this->mem->advectee().copy()(ii, jj);
      gf.resize(ii, jj);
      gf = this->mem->g_factor().copy()(ii, jj);

      //output stats
      ofs << "emin = " << (min(fc) - min(ic)) / max(ic)               << std::endl;
      ofs << "emax = " << (max(fc) - max(ic)) / max(ic)               << std::endl;
      ofs << "err0 = " << sqrt(sum(pow2(fc - ic) * gf)) / max(ic)     << std::endl;
      ofs << "err1 = " << sum(gf * fc) / sum(gf * ic) - 1             << std::endl;
      ofs << "err2 = " << sum(gf * pow2(fc)) / sum(gf * pow2(ic)) - 1 << std::endl;
    }
  }
};

