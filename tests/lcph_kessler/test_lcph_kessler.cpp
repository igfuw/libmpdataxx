#include <advoocat/solvers/mpdata_2d.hpp>
#include <advoocat/bcond/cyclic_2d.hpp>
#include <advoocat/concurr/threads.hpp>
#include <advoocat/output/gnuplot.hpp>

#include <libcloudph++/bulk_kessler.hpp>

enum {rhod_th_ix, rhod_rv_ix}; // TODO: better non-global, within main()

using namespace advoocat;

template <typename real_t, int n_iters, int n_eqs = 2>
class icicle : public solvers::mpdata_2d<real_t, n_iters, n_eqs>
{
  using parent_t = solvers::mpdata_2d<real_t, n_iters, n_eqs>;
  real_t acth;

  // deals with initials supersaturation
  void hook_ante_loop()
  {
  }

  // 
  void hook_post_step()
  {
    parent_t::hook_post_step();

    auto rhod_rv   = this->state(rhod_rv_ix);
    auto rhod_th   = this->state(rhod_th_ix);

    rng_t &i = this->i, &j = this->j;
    
    rhod_rv(i,j) += .1; // TODO: call something from bulk_kessler.hpp
  }

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t acth = .001; // TODO: document it! 
  };

  // ctor
  icicle( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p), 
    acth(p.acth)
  {}  
};


template <class T>
void setopts(T &p, int nt, int n_iters)
{
  //p.outfreq = nt/10; // TODO
  p.gnuplot_with = "lines";
  p.gnuplot_zrange = p.gnuplot_cbrange = "[.5:2.5]";
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << n_iters << "_%s_%d.svg";
    p.gnuplot_output = tmp.str();
  }
  p.outvars = {{rhod_rv_ix, {.name = "\\rho_d r_v", .unit = "TODO"}}};
}

int main()
{
  int nx = 32, nz = 32, nt = 200;

  const int n_iters = 2;
  using real_t = double;
  using solver_t = output::gnuplot<icicle<real_t, n_iters>>;
  solver_t::params_t p;
  setopts(p, nt, n_iters);
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(nx, nz, p);
  //setup(slv);
  slv.advance(nt);
}
