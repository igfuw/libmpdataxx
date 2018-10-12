/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include <boost/math/constants/constants.hpp>
#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>

using namespace libmpdataxx;
using T = double;

const T pi = boost::math::constants::pi<T>();

template <int opts_arg, int opts_iters, int sdiff_arg>
void test(int np, const std::string &name)
{
  std::string test_name = name + "_" + std::to_string(np);
  std::cout << "Calculating: " << test_name << std::endl;

  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { opts = opts_arg };
    enum { prs_order = 2 };
    enum { n_dims = 2 };
    enum { n_eqns = 2 };
    enum { rhs_scheme = solvers::trapez };
    enum { prs_scheme = solvers::cr };
    enum { sgs_scheme = solvers::dns };
    enum { stress_diff = sdiff_arg };
    struct ix { enum {
      u, v,
      vip_i=u, vip_j=v, vip_den=-1
    }; };
  
    enum { hint_norhs = opts::bit(ix::u) | opts::bit(ix::v)}; 
  }; 

  using ix = typename ct_params_t::ix;

  using solver_t = output::hdf5_xdmf<solvers::mpdata_rhs_vip_prs_sgs<ct_params_t>>;
  typename solver_t::rt_params_t p;

  p.di = 1.0 / (np - 1);
  p.dj = 1.0 / (np - 1);
  p.n_iters = opts_iters;
  p.prs_tol = 1e-8;
  p.grid_size = {np, np};
  p.eta = 1.0 / 20000;
  
  T time = 2.0;
  p.dt = 0.25 * p.di;
  int nt = time / p.dt;
  time = nt * p.dt;

  p.outfreq = nt / 8;
  p.outdir = "out_" + test_name;
  p.outvars = {
    {ix::u,   {"u",   "m/s"}}, 
    {ix::v,   {"v",   "m/s"}}
  };

  libmpdataxx::concurr::threads<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > slv(p);

  T dlta = 100;
  {
    blitz::firstIndex i;
    blitz::secondIndex j;

    slv.advectee(ix::u) = where(j * p.dj <= 0.5, tanh((j * p.dj - 0.25) * dlta), tanh((0.75 - j * p.dj) * dlta));
    slv.advectee(ix::v) = 0.05 * sin(2 * pi * i * p.di); 
  }

  slv.advance(nt); 
}

int main()
{
  std::vector<int> nps = {129};

  for (auto np : nps)
  {
    //test<opts::iga | opts::fct, 2, solvers::compact>(np, "iga_fct");
    //test<opts::iga | opts::tot | opts::fct, 2, solvers::compact>(np, "iga_tot_fct");
    test<opts::iga | opts::div_2nd | opts::div_3rd | opts::fct, 2, solvers::compact>(np, "iga_div3_fct");
  }
}
