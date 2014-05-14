//<listing-1>
#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

int main()
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqns = 1 };
  };

  // solver choice
  using slv_t = solvers::mpdata<ct_params_t>;

  // output choice
  using slv_out_t = output::gnuplot<slv_t>;
  
  // concurency choice
  using run_t = concurr::serial<
    slv_out_t, bcond::open, bcond::open
  >;         //left bcond   //right bcond

  // run-time parameters
  typename slv_out_t::rt_params_t p;

  int nx = 101, nt = 100;
  ct_params_t::real_t dx = 0.1;

  p.grid_size = { nx };
  p.outfreq = 20; 
 
  // instantiation
  run_t run(p);

  // initial condition
  blitz::firstIndex i;
  // Witch of Agnesi with a=.5 
  run.advectee() = -.5 + 1 / (
    pow(dx*(i - (nx-1)/2.), 2) + 1
  );
  // Courant number
  run.advector() = .5;

  // integration
  run.advance(nt);
}
//</listing-1>
