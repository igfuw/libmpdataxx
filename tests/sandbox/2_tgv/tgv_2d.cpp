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

struct err_t
{
  T Li, L2;
};

template <int opts_arg, int opts_iters, int sdiff_arg>
err_t test(int np)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { opts = opts_arg };
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
  
  using solver_t = libmpdataxx::solvers::mpdata_rhs_vip_prs_sgs<ct_params_t>;
  typename solver_t::rt_params_t p;

  p.di = 2 * pi / (np - 1);
  p.dj = 2 * pi / (np - 1);
  p.dt = 0.05 * p.di;
  p.n_iters = opts_iters;
  p.prs_tol = 1e-9;
  p.grid_size = {np, np};
  p.eta = 0.2;

  T time = 0.25 * pi;
  int nt = time / p.dt;
  time = nt * p.dt;

  libmpdataxx::concurr::threads<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > slv(p);

  decltype(slv.advectee(ix::u)) exact_u(slv.advectee(ix::u).shape());
  {
    blitz::firstIndex i;
    blitz::secondIndex j;

    slv.advectee(ix::u) =    cos(p.di * i) * sin(p.dj * j);
    slv.advectee(ix::v) =   -sin(p.di * i) * cos(p.dj * j); 

    exact_u =  cos(p.di * i) * sin(p.dj * j) * exp(-2 * p.eta * time);
  }

  slv.advance(nt); 

  err_t res;
  res.L2 = sqrt(sum(pow2(slv.advectee(ix::u) - exact_u))) / sqrt(sum(pow2(exact_u)));
  res.Li = max(abs(slv.advectee(ix::u) - exact_u)) / max(abs(exact_u));

  return res;
};

template <int opts_arg, int opts_iters>
bool check()
{
    // normal differencing of diffusive laplacian
    err_t err_coarse_nrml = test<opts_arg, opts_iters, solvers::normal>(33);
    err_t err_fine_nrml = test<opts_arg, opts_iters, solvers::normal>(65);

    int ord_Li_nrml = std::round(std::log2(err_coarse_nrml.Li / err_fine_nrml.Li));
    int ord_L2_nrml = std::round(std::log2(err_coarse_nrml.L2 / err_fine_nrml.L2));

    // compact differencing of diffusive laplacian
    err_t err_coarse_cmpct = test<opts_arg, opts_iters, solvers::compact>(33);
    err_t err_fine_cmpct = test<opts_arg, opts_iters, solvers::compact>(65);

    int ord_Li_cmpct = std::round(std::log2(err_coarse_cmpct.Li / err_fine_cmpct.Li));
    int ord_L2_cmpct = std::round(std::log2(err_coarse_cmpct.L2 / err_fine_cmpct.L2));

    return ord_Li_nrml  == 2 && ord_L2_nrml  == 2 &&
           ord_Li_cmpct == 2 && ord_L2_cmpct == 2 &&
           err_fine_cmpct.L2 < err_fine_nrml.L2   && err_fine_cmpct.Li < err_fine_nrml.Li;
}

int main()
{
  if (!check<opts::iga, 2>()) throw std::runtime_error("iga");
}
