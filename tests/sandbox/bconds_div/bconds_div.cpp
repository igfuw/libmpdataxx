// test to check if we obtain non-divergent advector field with requested precision
// using different sets of boundary conditions, setup based on the boussinesq test
// from the paper suite
#include <libmpdata++/solvers/boussinesq.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include "bconds_div.hpp"

using namespace libmpdataxx;

template <bcond::bcond_e bcond_h, bcond::bcond_e bcond_v>
void test(const std::string &error_str)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqns = 3 };
    enum { rhs_scheme = solvers::trapez };
    enum { prs_scheme = solvers::cr };
    struct ix { enum {
      u, w, tht, 
      vip_i=u, vip_j=w, vip_den=-1
    }; };
  }; 
  
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  
  const int r0 = 250; 
  const int nx = 201, ny = 201, nt = 10;
  typename ct_params_t::real_t Tht_ref = 300;

  using slv_t = bconds_div<ct_params_t>;

  typename slv_t::rt_params_t p;
  p.dt = 7.5;
  p.di = p.dj = 10.; 
  p.Tht_ref = Tht_ref; 
  p.prs_tol = 1e-10;
  p.grid_size = {nx, ny};
  p.error_str = error_str;

  libmpdataxx::concurr::threads<
    slv_t, 
    bcond_h, bcond_h,
    bcond_v, bcond_v
  > slv(p);

  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    
    slv.sclr_array("tht_e") = Tht_ref;
    slv.advectee(ix::tht) = Tht_ref + where(
      // if
      pow(i * p.di - 4    * r0 , 2) + 
      pow(j * p.dj - 1.04 * r0 , 2) <= pow(r0, 2), 
      // then
      .5, 
      // else
      0
    );
    slv.advectee(ix::u) = 0; 
    slv.advectee(ix::w) = 0; 
  }

  slv.advance(nt);  
}

int main() 
{
  test<bcond::cyclic, bcond::cyclic>("cyclic_cyclic");
  test<bcond::open  , bcond::cyclic>("open_cyclic");
  test<bcond::open  , bcond::rigid >("open_rigid");
  test<bcond::cyclic, bcond::rigid >("cyclic_rigid");
};
