/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/concurr/serial.hpp>

using namespace libmpdataxx;

template<int vip_vab_v>
void test(double coeff, int nt, double tol, const std::string& name)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqns = 2 };
    enum { rhs_scheme = solvers::euler_a };
    enum { vip_vab = vip_vab_v};
    enum { prs_scheme = solvers::cr};
    struct ix { enum {
      u, w,
      vip_i=u, vip_j=w, vip_den=-1
    }; };
    enum { hint_norhs = opts::bit(ix::u) | opts::bit(ix::w)}; 
  }; 

  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;

  using slv_t = solvers::mpdata_rhs_vip_prs<ct_params_t>;
  typename slv_t::rt_params_t p;
  p.dt = 1;
  p.di = p.dj = 1; 

  p.prs_tol = 1e-6;
  p.grid_size = {4, 4};

  libmpdataxx::concurr::serial<
    slv_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > slv(p);

  real_t u0 = 0.2;
  real_t w0 = 0.2;

  slv.advectee(ix::u) = u0; 
  slv.advectee(ix::w) = w0; 
    
  slv.vab_coefficient() = coeff;
  
  real_t ur = 0.15;
  real_t wr = 0.25;
  slv.vab_relaxed_state(0) = ur;
  slv.vab_relaxed_state(1) = wr;

  slv.advance(nt);
  real_t t = nt * p.dt;

  auto diff_u = max(abs(slv.advectee(ix::u) - ur));
  auto err_u = max(abs(slv.advectee(ix::u) - (ur + (u0 - ur) * exp(-coeff * t))));

  auto diff_w = max(abs(slv.advectee(ix::w) - wr));
  auto err_w = max(abs(slv.advectee(ix::w) - (wr + (w0 - wr) * exp(-coeff * t))));
 
  if (nt == 1)
  {
    if (diff_u > 0)
    { 
      std::cerr << diff_u << std::endl;
      throw(std::runtime_error("u " + name));
    }
    if (diff_w > 0)
    {
      std::cerr << diff_w << std::endl;
      throw(std::runtime_error("w " + name));
    }
  }
  else
  {
    if (err_u > tol)
    {
      std::cerr << err_u << std::endl;
      throw(std::runtime_error("u " + name));
    }
    if (err_w > tol)
    {
      std::cerr << err_w << std::endl;
      throw(std::runtime_error("w " + name));
    }
  }
}

int main() 
{
#if defined(USE_MPI)
  // we will instantiate many solvers, so we have to init mpi manually, 
  // because solvers will not know should they finalize mpi upon destruction
  MPI::Init_thread(MPI_THREAD_SERIALIZED);
#endif
  test<solvers::expl>(1, 1, 0, "explA");
  test<solvers::expl>(0.01, 100, 1e-4, "explB");
  test<solvers::impl>(2, 1, 0, "implA");
  test<solvers::impl>(0.01, 100, 2e-7, "imblB");
#if defined(USE_MPI)
  MPI::Finalize();
#endif
};
