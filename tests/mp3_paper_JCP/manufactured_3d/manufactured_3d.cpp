#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include "../../sandbox/convergence_spacetime/convergence_3d_spacetime.hpp" // reusing convergence_3d_spacetime test

using namespace libmpdataxx;

template <int opts_arg, int opts_iters>
void test(const int np, const std::string &name)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 3 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
    enum { var_dt = true };
    enum { out_intrp_ord = 2 };
  };

  typename ct_params_t::real_t
    pi = boost::math::constants::pi<typename ct_params_t::real_t>(),
    time = 1.0,
    dx = 2.0 * pi / (np - 1),
    dy = 2.0 * pi / (np - 1),  
    dz = 2.0 * pi / (np - 1), 
    dt = 0.5 * dx;

  using slv_t = output::hdf5_xdmf<convergence_3d_spacetime<ct_params_t>>;

  typename slv_t::rt_params_t p;
  p.n_iters = opts_iters; 
  p.grid_size = {np, np, np};
  p.di = dx;
  p.dj = dy;
  p.dk = dz;
  p.dt = dt;
  p.max_courant = 0.5;

  p.outfreq = 1.0;
  p.outdir = "out_" + name + "_" + std::to_string(np);

  concurr::threads<
    slv_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p); 

  run.advectee() = 8;

  // just to fool the divergence check
  run.advector(0) = 0;
  run.advector(1) = 0;
  run.advector(2) = 0;
  
  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;

  run.g_factor() = exp(cos(i * dx) + cos(j * dy) + cos(k * dz));

  run.advance(time);
}

int main()
{
  std::vector<int> nps = {9, 17, 33, 65, 129};

  for (const auto np : nps)
  {
    {
      enum { opts = opts::nug | opts::tot | opts::dfl};
      enum { opts_iters = 3};
      test<opts, opts_iters>(np, "Mp3cc");
    }
    
    {
      enum { opts = opts::nug | opts::div_2nd | opts::div_3rd};
      enum { opts_iters = 2};
      test<opts, opts_iters>(np, "Mp3");
    }
  }
}
