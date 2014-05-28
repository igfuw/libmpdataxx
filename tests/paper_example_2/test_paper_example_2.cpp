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
//<listing-1>
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
//</listing-1>
    enum { opts = opts_arg };
  };

  using sim_t = output::gnuplot<
    solvers::mpdata<ct_params_t>
  >;
  typename sim_t::rt_params_t p;

//<listing-2>
  int nx = 601, nt = 1200;
  // run-time parameters
  p.grid_size = { nx };
  p.outfreq = nt; 
  p.outvars = {
    {0, {.name = "", .unit = "1"}},
    {1, {.name = "", .unit = "1"}}
  };
//</listing-2>
  p.gnuplot_output = filename; 
  p.gnuplot_command = "plot";
  p.gnuplot_with = "histeps";
  p.gnuplot_yrange = "[-1.25:4.25]";
  p.gnuplot_fontsize = "17";

  // instantiation
  concurr::serial<sim_t, bcond::cyclic, bcond::cyclic> run(p);

//<listing-3>
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
//</listing-3>

  // integration
  run.advance(nt);
}

int main()
{
  {
//<listing-4>
    enum { opts = opts::abs };
//</listing-4>
    test<opts>("out_abs.svg");
  }
  {
//<listing-5>
    enum { opts = opts::iga };
//</listing-5>
    test<opts>("out_iga.svg");
  }
  {
//<listing-6>
    enum { opts = opts::iga | opts::tot };
//</listing-6>
    test<opts>("out_iga_tot.svg");
  }
  {
//<listing-7>
    enum { opts = opts::iga | opts::fct };
//</listing-7>
    test<opts>("out_iga_fct.svg");
  }
  {
//<listing-8>
    enum { opts = opts::iga | opts::tot | opts::fct };
//</listing-8>
    test<opts>("out_iga_tot_fct.svg");
  }
  {
    enum { opts = opts::iga | opts::fct | opts::khn };
    test<opts>("out_iga_fct_khn.svg");
  }

}
