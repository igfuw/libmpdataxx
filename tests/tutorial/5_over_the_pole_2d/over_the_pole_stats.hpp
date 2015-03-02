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

  real_t last_timestep;

  blitz::Range ii, jj;
 
  std::ofstream ofs;
  
  template <class T>
  void write_stat(const std::string& stat_info,
                  const T stat_val,
                  const bool protect_from_neg_zero = true,
                  const int precision = 8)
  {
      const auto neg_zero = "-0." + std::string(precision, '0');

      std::stringstream ss;
      ss << std::fixed << std::setprecision(precision) << stat_val;
      
      std::string stat_val_s;
      if (protect_from_neg_zero)
        stat_val_s = (ss.str() == neg_zero ? ss.str().erase(0, 1) : ss.str());
      else
        stat_val_s = ss.str();

      ofs << stat_info << stat_val_s << std::endl;
  }

  void hook_ante_loop(const int nt)
  {
    parent_t::hook_ante_loop(nt);
    if (this->rank != 0) return;

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
    ic = this->mem->advectee()(ii, jj);

    // at the end of the test true solution is equal to the inital state
    ofs << std::endl;
    write_stat("timestep      = ", 0);
    write_stat("min(solution) = ", min(ic), false);
    write_stat("max(solution) = ", max(ic));
    ofs << " " << std::endl;
  }

  void hook_post_step()
  {
    parent_t::hook_post_step();
    this->mem->barrier();
    if (this->rank != 0) return;
    if (this->timestep == last_timestep) 
    { 
      // final condition
      fc.resize(ii, jj);
      fc = this->mem->advectee()(ii, jj);
      gf.resize(ii, jj);
      gf = this->mem->g_factor()(ii, jj);

      //output stats
      write_stat("timestep = ", this->timestep);
      write_stat("emin = ", (min(fc) - min(ic)) / max(ic));
      write_stat("emax = ", (max(fc) - max(ic)) / max(ic));
      write_stat("err0 = ", sqrt(sum(pow2(fc - ic) * gf)) / max(ic));
      write_stat("err1 = ", sum(gf * fc) / sum(gf * ic) - 1);
      write_stat("err2 = ", sum(gf * pow2(fc)) / sum(gf * pow2(ic)) - 1);
    }
  }
};
