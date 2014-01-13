//<listing>
#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

int main()
{
  // compile-time parameters
  struct ct_params_t
  {
    using real_t = float;
    enum { n_dims = 1 };
    enum { n_eqs = 1 };
    enum { opts = 0 };
  };

  // solver choice
  using slv_t = solvers::mpdata<ct_params_t>;

  // output choice
  using sim_t = output::gnuplot<slv_t>;
  
  // concurency choice
  using run_t = concurr::serial<
    sim_t, bcond::cyclic
  >;

  // run-time parameters
  typename sim_t::rt_params_t p;

  int nx = 128, nt = 128;
  ct_params_t::real_t dx = 0.128;

  p.span = { nx };
  p.n_iters = 2;
  p.outfreq = nt / 10; 
 
  // instantiation
  run_t run(p);

  // initial condition
  blitz::firstIndex i;
  // Witch of Agnesi with a=.5
  run.advectee() = 1 / (
    pow(dx*(.5+i - nx/2.), 2) + 1
  );
  // Courant number
  run.advector() = .5;

  // integration
  run.advance(nt);
}
//</listing>
