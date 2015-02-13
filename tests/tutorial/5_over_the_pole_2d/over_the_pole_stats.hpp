/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * errors for over the pole test
 */
#pragma once

template<class parent_t>
struct stats : public parent_t
{
  using parent_t::parent_t;
  using real_t = typename parent_t::real_t;

  typename parent_t::arr_t ic, fc, gf;

  real_t last_timestep, err0, err1, err2;

  blitz::Range ii, jj;
 
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
                                    //nlon
    ii = blitz::Range(0, this->mem->grid_size[0] - 2);  // due to boundary condition type not all points should be evaluated
    jj = blitz::Range(0, this->mem->grid_size[1] - 1);  // when calculating errors
                                    //nlat

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
