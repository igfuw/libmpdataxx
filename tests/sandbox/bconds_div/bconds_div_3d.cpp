// test to check if we obtain non-divergent advector field with requested precision
// using different sets of boundary conditions, setup based on the boussinesq test
// from the paper suite
#include <libmpdata++/solvers/boussinesq.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include "bconds_div.hpp"

using namespace libmpdataxx;

template <bcond::bcond_e bcond_x, bcond::bcond_e bcond_y, bcond::bcond_e bcond_z>
void test(const std::string &error_str)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 3 };
    enum { n_eqns = 4 };
    enum { rhs_scheme = solvers::trapez };
    enum { prs_scheme = solvers::cr };
    struct ix { enum {
      u, v, w, tht, 
      vip_i=u, vip_j=v, vip_k=w, vip_den=-1
    }; };
  }; 
  
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  
  const double r0 = 37.5; 
  const int nx = 31, ny = 31, nz = 31, nt = 10;
  typename ct_params_t::real_t Tht_ref = 300;

  using slv_t = bconds_div<ct_params_t>;

  typename slv_t::rt_params_t p;
  p.dt = 15;
  p.di = p.dj = p.dk = 10.; 
  p.Tht_ref = Tht_ref; 
  p.prs_tol = 1e-10;
  p.grid_size = {nx, ny, nz};
  p.error_str = error_str;

  libmpdataxx::concurr::threads<
    slv_t, 
    bcond_x, bcond_x,
    bcond_y, bcond_y,
    bcond_z, bcond_z
  > slv(p);

  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;
    
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
    slv.advectee(ix::v) = 0; 
    slv.advectee(ix::w) = 0; 
  }

  slv.advance(nt);  
}

int main() 
{
  test<bcond::cyclic, bcond::cyclic, bcond::cyclic>("cyclic_cyclic_cyclic");
  test<bcond::cyclic, bcond::cyclic, bcond::rigid >("cyclic_cyclic_rigid" );
  test<bcond::open  , bcond::open  , bcond::rigid >("open_open_rigid"     );
  test<bcond::open  , bcond::rigid , bcond::cyclic>("open_rigid_cyclic"   );
};
