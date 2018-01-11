/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

#include "../mp3_paper_2018_JCP/reversing_deform/test_def.hpp"

using namespace libmpdataxx;

int main()
{
  const bool var_dt = true;
  const T max_cfl = 0.90;

  std::vector<int> nys = {120, 240};
  {
    enum { opts = opts::nug};
    const int opts_iters = 2;
    for (const auto ny : nys) test<var_dt, opts, opts_iters>("nug_i2", ny, max_cfl);
  }
  
  {
    enum { opts = opts::nug | opts::iga | opts::fct};
    const int opts_iters = 2;
    for (const auto ny : nys) test<var_dt, opts, opts_iters>("nug_iga_fct_i2", ny, max_cfl);
  }
}
