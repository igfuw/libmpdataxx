/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/solvers/shallow_water.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/hdf5.hpp>
using namespace libmpdataxx; 

using real_t = double;

struct intcond
{
  real_t operator()(const real_t &x) const
  {
    return 
      std::abs(x) <= 1 // if
      ? 1 - x*x        // then
      : 0;             // else
  }
  BZ_DECLARE_FUNCTOR(intcond);
};

template<int opts_arg>
void test(const std::string& outdir) 
{
  // compile-time parameters
  // enum { hint_noneg = opts::bit(ix::h) };  // TODO: reconsider?
  struct ct_params_t : ct_params_default_t
  {
    using real_t = ::real_t;
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
  
    // options
    enum { opts = opts_arg | opts::dfl };
    enum { rhs_scheme = solvers::trapez };
  
    // indices
    struct ix { 
      enum { qx, h };
      enum { vip_i=qx, vip_den=h };
    }; 
  
    // hints
    enum { hint_norhs = opts::bit(ix::h) }; 
  };

  const int 
    nt = 300,
    outfreq = 100;

  using ix = typename ct_params_t::ix;

  // solver choice
  using slv_out_t =
    output::hdf5<
      solvers::shallow_water<ct_params_t>
    >;
 
  // run-time parameters
  typename slv_out_t::rt_params_t p;
  p.dt = .01;
  p.di = .05;
  p.grid_size = { int(16 / p.di) };
  p.g = 1;
  p.vip_eps = 1e-8; 
  p.outvars[ix::qx].name = "qx";
  p.outvars[ix::h].name  = "h";
  p.outdir = outdir;
  p.outfreq = outfreq; 
 
  // instantiation
  concurr::serial<
    slv_out_t, 
    bcond::open, bcond::open
  > run(p);

  // initial condition
  {
    blitz::firstIndex i;
    run.advectee(ix::h) = intcond()(p.di * i - (p.grid_size[0]-1) * p.di / 2);
  }
  run.advectee(ix::qx) = 0;

  // integration
  run.advance(nt); 
}

int main()
{
#if defined(USE_MPI)
  // we will instantiate many solvers, so we have to init mpi manually, 
  // because solvers will not know should they finalize mpi upon destruction
  MPI::Init_thread(MPI_THREAD_MULTIPLE);
#endif
  test<opts::abs | opts::fct >("1d_fct_abs");
  test<opts::iga | opts::fct >("1d_fct_iga");
#if defined(USE_MPI)
  MPI::Finalize();
#endif
}
