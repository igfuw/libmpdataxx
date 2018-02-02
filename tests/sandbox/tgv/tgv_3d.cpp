/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include <boost/math/constants/constants.hpp>

using namespace libmpdataxx;
using boost::math::constants::pi;

void test(const std::string &dirname, double rey)
{
  const int nx = 65, ny = 65, nz = 65, nt = 2000;

  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 3 };
    enum { opts = opts::fct | opts::iga};
    enum { n_eqns = 3 };
    enum { rhs_scheme = solvers::trapez };
    enum { sgs_scheme = solvers::dns };
    enum { stress_diff = solvers::pade };
    enum { prs_scheme = solvers::cr };
    struct ix { enum {
      u, v, w,
      vip_i=u, vip_j=v, vip_k=w, vip_den=-1
    }; };
  
    enum { hint_norhs = opts::bit(ix::u) | opts::bit(ix::v) | opts::bit(ix::w)}; 
  }; 

  using ix = typename ct_params_t::ix;

  using solver_t = output::hdf5_xdmf<libmpdataxx::solvers::mpdata_rhs_vip_prs_sgs<ct_params_t>>;

  solver_t::rt_params_t p;

  p.n_iters = 2;

  p.eta = 1.0 / rey;

  p.dt = 0.005;
  p.di = 2 * pi<ct_params_t::real_t>() / (nx - 1);
  p.dj = 2 * pi<ct_params_t::real_t>() / (ny - 1);
  p.dk = 2 * pi<ct_params_t::real_t>() / (nz - 1);

  p.outfreq = 50; 
  p.outwindow = 3;
  p.outvars = {
    {ix::u,   { "u",   "m/s"}}, 
    {ix::v,   { "v",   "m/s"}}, 
    {ix::w,   { "w",   "m/s"}}, 
  };
  p.outdir = dirname;
  p.prs_tol = 1e-7;
  p.grid_size = {nx, ny, nz};

  libmpdataxx::concurr::threads<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > slv(p);

  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    slv.advectee(ix::u) =  sin(p.di * i) * cos(p.dj * j) * cos(p.dk * k);
    slv.advectee(ix::v) = -cos(p.di * i) * sin(p.dj * j) * cos(p.dk * k); 
    slv.advectee(ix::w) = 0.0; 
  }

  slv.advance(nt); 
};

int main()
{
  test("rey=800", 800);
  test("rey=100", 100);
}
