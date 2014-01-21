#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

int main()
{
//<listing-1>
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqs = 2 };
    enum { opts = formulae::opts::iga };
  };
//</listing-1>

  using sim_t = output::gnuplot<
    solvers::mpdata<ct_params_t>
  >;
  typename sim_t::rt_params_t p;

  int nx = 500, nt = 1600;

//<listing-2>
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
//</listing-2>

  // instantiation
  concurr::serial<sim_t, bcond::cyclic> run(p);

//<listing-3>
  // initial condition
  blitz::firstIndex i;
  run.advectee(0) = where(
    i <= 75 || i >= 125, // if
    2,                   // then
    4                    // else
  ); 
  run.advectee(1) = where(
    i <= 75 || i >= 125, // if 
    -1,                  // then
    1                    // else
  ); 
  run.advector() = -.5;  // Courant
//</listing-3>

  // integration
  run.advance(nt);
}
