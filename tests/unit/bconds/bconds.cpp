/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <random>

using namespace libmpdataxx;

int main() 
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 3 };
    enum { n_eqns = 3 };
    enum { rhs_scheme = solvers::trapez };
    enum { prs_scheme = solvers::cr};
    struct ix { enum {
      u, v, w,
      vip_i=u, vip_j=v, vip_k=w, vip_den=-1
    }; };
    enum { hint_norhs = opts::bit(ix::u) | opts::bit(ix::v) | opts::bit(ix::w)}; 
  }; 

  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;

  const int nx = 17, ny = 9, nz = 13, nt = 100;

  using slv_t = solvers::mpdata_rhs_vip_prs<ct_params_t>;

  slv_t::rt_params_t p;

  p.di = 1;
  p.dj = 1;
  p.dk = 1;
  p.dt = 0.1;
  p.prs_tol = 1e-10;
  p.grid_size = {nx, ny, nz};

  libmpdataxx::concurr::threads<
    slv_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::open, bcond::open,
    bcond::rigid, bcond::rigid
  > slv(p);

  slv.advectee(ix::u) = 1;
  slv.advectee(ix::w) = 0;

  blitz::thirdIndex k;
  slv.advectee(ix::v) = k / real_t(nz - 1); 
  
  // save initial v profile
  auto all = blitz::Range::all();
  auto vprof = slv.advectee_global(ix::v)(all, 0, all);

  std::mt19937 gen(1);
  std::uniform_real_distribution<> dis(-1, 1);
  real_t amp = 1e-6;
  auto rand = std::bind(dis, gen);

  // add random perturbation in the interior
  std::array<decltype(ix::u), 3> velo = {ix::u, ix::v, ix::w};
  std::for_each(velo.begin(), velo.end(), [&](decltype(ix::u) vel)
  {
    decltype(slv.advectee(vel)) prtrb(slv.advectee_global(vel).shape()); // array to store perturbation
    std::generate(prtrb.begin(), prtrb.end(), [&] () {return amp * rand();}); // fill it, TODO: is it officialy stl compatible?
    // no perturbation at the edges
    prtrb(0, all, all) = 0;
    prtrb(all, 0, all) = 0;
    prtrb(all, all, 0) = 0;
    prtrb(nx-1, all, all) = 0;
    prtrb(all, ny-1, all) = 0;
    prtrb(all, all, nz-1) = 0;
    // apply the pperturbation
    prtrb += slv.advectee_global(vel);
    slv.advectee_global_set(prtrb, vel);
  } );

  slv.advance(nt);

  auto err_cyclic = max(abs(slv.advectee_global(ix::u)(0, all, all) - slv.advectee_global(ix::u)(nx - 1, all, all)));
  // TODO: why MPI run has larger error?
  const real_t cyclic_tolerance = 
  #if defined(USE_MPI)
    1e-6;
  #else
    1e-10;
  #endif
  if (err_cyclic > cyclic_tolerance)
  {
    std::cerr << err_cyclic << std::endl;
    throw std::runtime_error("cyclic");
  }

  auto err_rigid1 = max(abs(slv.advectee_global(ix::w)(all, all, 0)));
  if (err_rigid1 > 0)
  {
    std::cerr << err_rigid1 << std::endl;
    throw std::runtime_error("rigid1");
  }
  auto err_rigid2 = max(abs(slv.advectee_global(ix::w)(all, all, nz - 1)));
  if (err_rigid2 > 0)
  {
    std::cerr << err_rigid2 << std::endl;
    throw std::runtime_error("rigid2");
  }

  auto err_open1 = max(abs(slv.advectee_global(ix::v)(all, 0, all) - vprof));
  if (err_open1 > 0)
  {
    std::cerr << err_open1 << std::endl;
    throw std::runtime_error("open1");
  }
  auto err_open2 = max(abs(slv.advectee_global(ix::v)(all, ny - 1, all) - vprof));
  if (err_open2 > 0)
  {
    std::cerr << err_open2 << std::endl;
    throw std::runtime_error("open2");
  }
}
