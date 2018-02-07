#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/gnuplot.hpp>

using namespace libmpdataxx;

template <int opts_arg>
void test(const std::string filename)
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
    enum { opts = opts_arg };
  };

  using sim_t = output::gnuplot<
    solvers::mpdata<ct_params_t>
  >;
  typename sim_t::rt_params_t p;

  int nx = 601, nt = 1200;
  // run-time parameters
  p.grid_size = { nx };
  p.outfreq = nt; 
  p.outvars = {
    {0, {"", "1"}},
    {1, {"", "1"}}
  };
  p.gnuplot_output = filename; 
  p.gnuplot_command = "plot";
  p.gnuplot_with = "histeps";
  p.gnuplot_yrange = "[-1.25:4.25]";
  p.gnuplot_fontsize = "17";

  // instantiation
  concurr::serial<sim_t, bcond::cyclic, bcond::cyclic> run(p);

  // initial condition
  blitz::firstIndex i;
  run.advectee(0) = where(
    i <= 75 || i >= 125,   // if
    2,                     // then
    4                      // else
  ); 
  run.advectee(1) = where(
    i <= 75 || i >= 125,   // if 
    -1,                    // then
    1                      // else
  ); 
  run.advector() = -.75;  // Courant

  // integration
  run.advance(nt);
}

int main()
{
  {
    enum { opts = opts::abs };
    test<opts>("out_abs.svg");
  }
  {
    enum { opts = opts::iga };
    test<opts>("out_iga.svg");
  }
  {
    enum { opts = opts::iga | opts::tot };
    test<opts>("out_iga_tot.svg");
  }
  {
    enum { opts = opts::iga | opts::fct };
    test<opts>("out_iga_fct.svg");
  }
  {
    enum { opts = opts::iga | opts::tot | opts::fct };
    test<opts>("out_iga_tot_fct.svg");
  }
}
