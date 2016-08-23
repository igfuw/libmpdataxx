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
  auto vprof = slv.advectee(ix::v)(all, 0, all);

  std::mt19937 gen(1);
  std::uniform_real_distribution<> dis(-1, 1);
  real_t amp = 1e-6;
  
  // add random perturbation in the interior
  for (int i = 1; i < nx - 1; ++i)
  {
    for (int j = 1; j < ny - 1; ++j)
    {
      for (int k = 1; k < nz - 1; ++k)
      {
        slv.advectee(ix::u)(i, j, k) += amp * dis(gen);
        slv.advectee(ix::v)(i, j, k) += amp * dis(gen);
        slv.advectee(ix::w)(i, j, k) += amp * dis(gen);
      }
    }
  }

  slv.advance(nt);

  auto err_cyclic = max(abs(slv.advectee(ix::u)(0, all, all) - slv.advectee(ix::u)(nx - 1, all, all)));
  if (err_cyclic > 1e-10)
  {
    std::cerr << err_cyclic << std::endl;
    throw std::runtime_error("cyclic");
  }

  auto err_rigid1 = max(abs(slv.advectee(ix::w)(all, all, 0)));
  if (err_rigid1 > 0)
  {
    std::cerr << err_rigid1 << std::endl;
    throw std::runtime_error("rigid1");
  }
  auto err_rigid2 = max(abs(slv.advectee(ix::w)(all, all, nz - 1)));
  if (err_rigid2 > 0)
  {
    std::cerr << err_rigid2 << std::endl;
    throw std::runtime_error("rigid2");
  }

  auto err_open1 = max(abs(slv.advectee(ix::v)(all, 0, all) - vprof));
  if (err_open1 > 0)
  {
    std::cerr << err_open1 << std::endl;
    throw std::runtime_error("open1");
  }
  auto err_open2 = max(abs(slv.advectee(ix::v)(all, ny - 1, all) - vprof));
  if (err_open2 > 0)
  {
    std::cerr << err_open2 << std::endl;
    throw std::runtime_error("open2");
  }
}
