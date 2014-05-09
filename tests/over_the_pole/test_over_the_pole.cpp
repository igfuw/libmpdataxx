/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <cmath>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
using namespace libmpdataxx;

template <int opts_arg, int opts_iters>
void test(const std::string filename)
{
  enum {x, y};
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqs = 1 };
    enum { opts = opts_arg };
  };
  
  int
    nlon = 129,
    nlat = 64,
    nt = 5120;
  
  typename ct_params_t::real_t 
    pi = boost::math::constants::pi<typename ct_params_t::real_t>(),
    dlmb = 2 * pi / (nlon - 1),
    dphi = pi / nlat;

  using slv_out_t = output::hdf5_xdmf<solvers::mpdata<ct_params_t>>;
  typename slv_out_t::rt_params_t p;

  p.n_iters = opts_iters; 
  p.span = {nlon, nlat};

  p.outfreq = nt; 
  p.outvars[0].name = "psi";
  p.outdir = filename;

//<listing-2>
  concurr::serial<
    slv_out_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::polar, bcond::polar
  > run(p); 
//</listing-2>
  
  typename ct_params_t::real_t
    r = 7 * dlmb,
    x0 = 3 * pi / 2,
    y0 = 0;

  blitz::firstIndex i;
  blitz::secondIndex j;

  decltype(run.advectee()) 
    tmp(run.advectee().extent());

  tmp = 2 * (  blitz::pow2(blitz::cos(dphi * (j + 0.5) - pi / 2) * blitz::sin((dlmb * i - x0) / 2))
             + blitz::pow2(blitz::sin((dphi * (j + 0.5) - pi / 2 - y0) / 2))                     );

  run.advectee() = where(
    tmp - pow(r, 2) <= 0,                  //if
    1 - sqrt(tmp) / r,   //then
    0.                                     //else
  );
  
  typename ct_params_t::real_t
    udt = 2 * pi / nt,
    b = -pi / 2;

  run.advector(x) = dlmb * udt * (cos(b) * blitz::cos((j+0.5) * dphi - pi / 2) + sin(b) * blitz::sin((j+0.5) * dphi - pi / 2) * blitz::cos((i+0.5) * dlmb));
  
  run.advector(y) = -dlmb * udt * sin(b) * blitz::sin((i) * dlmb)* blitz::cos((j+1) * dphi - pi / 2);
  
//<listing-3>
  run.g_factor() = dlmb * dphi *
    blitz::cos(dphi * (j + 0.5) - pi / 2);
//</listing-3>

  run.advance(nt);
}

int main()
{
  {
//<listing-1>
    enum { opts = formulae::opts::nug };
//</listing-1>

    enum { opts = formulae::opts::nug | formulae::opts::iga | formuale::opts::fct };
    const int opts_iters = 2;
    test<opts, opts_iters>("default");
    
    enum { opts = formulae::opts::nug | formulae::opts::tot | formuale::opts::fct };
    const int opts_iters = 3;
    test<opts, opts_iters>("default");
  }
}
