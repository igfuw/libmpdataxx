#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

int main()
{
  // compile-time parameters
  struct ct_params_t 
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqs = 2 };
    enum { opts = formulae::opts::iga };
  };

  using sim_t = output::gnuplot<
    solvers::mpdata<ct_params_t>
  >;
  typename sim_t::rt_params_t p;

  int nx = 500, nt = 1600;
  ct_params_t::real_t min[2] = {2, -1}, max[2] = {4, 1};

  // run-time parameters
  p.span = { nx };

  p.outfreq = nt; // diplays initial condition and the final state
  p.outvars = {
    {0, {.name = "single-sign signal", .unit = "1"}},
    {1, {.name = "variable-sign signal", .unit = "1"}}
  };
  p.gnuplot_command = "plot";
  p.gnuplot_with = "steps";
  p.gnuplot_yrange = "[-1.25:4.25]";

  // instantiation
  concurr::serial<sim_t, bcond::cyclic> run(p);

  // initial condition
  blitz::firstIndex i;
  int width = 50, center = 100;
  run.advectee(0) = where(i <= center-width/2 || i >= center+width/2, min[0], max[0]); 
  run.advectee(1) = where(i <= center-width/2 || i >= center+width/2, min[1], max[1]); 
  run.advector() = -.5; 

  // integration
  run.advance(nt);
}
