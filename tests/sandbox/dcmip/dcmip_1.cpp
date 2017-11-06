/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include "dcmip_1.hpp"

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

using namespace libmpdataxx;
using T = double;
const T pi = boost::math::constants::pi<T>();

template <int opts_arg, int opts_iters>
void test(const std::string filename, int ny)
{
  T R_d = 287.;
  T T_0 = 300.;
  T z_top = 12000.;
  T p_top = 254.944 * 100;
  T a = 6.37122e6;
  T g = 9.80616;
  T p_0 = 1000. * 100;
  T day = 86400.0;

  T H = R_d * T_0 / g;
  T rho_0 = p_0 / (R_d * T_0);

  T tau = 12 * day;
  T omg_0 = 23000 * pi / tau; // hPa / s 
  T b = 0.2;
  T x_c1 = 5 * pi / 6;
  T x_c2 = 7 * pi / 6;
  T y_ca = 0;
  T z_ca = 5000;
  T r_t =  a / 2;
  T z_t =  1000;

  const int nx = 2 * ny + 1;
  const int nz = ny / 3 + 1;

  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { var_dt = true };
    enum { n_dims = 3 };
    enum { n_eqns = 4 };
    enum { opts = opts_arg };
  };
  
  T dx = 2 * pi / (nx - 1), dy = pi / ny, dz = z_top / (nz - 1);
  
  T time = 12 * day;
  //T dt = dt_arg;
  //int nt = time / dt;

  using slv_out_t = 
      output::hdf5_xdmf<
        dcmip_1<ct_params_t>
    >;
  typename slv_out_t::rt_params_t p;

  p.n_iters = opts_iters; 
  p.grid_size = {nx, ny, nz};

  p.dt = 1;
  p.max_courant = 0.9;
  p.di = dx;
  p.dj = dy;
  p.dk = dz;

  p.outfreq = day; 
  p.outvars[0].name = "q1";
  p.outvars[1].name = "q2";
  p.outvars[2].name = "q3";
  p.outvars[3].name = "q4";
  p.outdir = filename;

  p.pi = pi;
  p.a = a;
  p.g = g;
  p.z_top = z_top;
  p.p_top = p_top;
  p.rho_0 = rho_0;
  p.p_0 = p_0;
  p.tau = tau;
  p.H = H;

  p.omg_0 = omg_0;
  p.b = b;

  concurr::threads<
    slv_out_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::polar, bcond::polar,
    bcond::rigid, bcond::rigid
  > run(p); 
  
  decltype(run.advectee()) solution(run.advectee().extent());

  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int k = 0; k < nz; ++k)
      {
        T x = i * dx;
        T y = (j + 0.5) * dy - pi / 2;
        T z = k * dz;

        T rho = rho_0 * exp(-z / H);

        auto r1 = a * acos(sin(y_ca) * sin(y) + cos(y_ca) * cos(y) * cos(x - x_c1));
        auto r2 = a * acos(sin(y_ca) * sin(y) + cos(y_ca) * cos(y) * cos(x - x_c2));
        
        auto d1 = std::min(1., pow(r1 / r_t, 2) + pow((z - z_ca) / z_t, 2));
        auto d2 = std::min(1., pow(r2 / r_t, 2) + pow((z - z_ca) / z_t, 2));
        
        run.advectee(0)(i, j, k) = 0.5 * (1 + cos(pi * d1)) + 0.5 * (1 + cos(pi * d2));

        run.advectee(1)(i, j, k) = 0.9 - 0.8 * pow(run.advectee(0)(i, j, k), 2);

        run.advectee(2)(i, j, k) = (d1 < 0.5 || d2 < 0.5) ? 1. : 0.;
        run.advectee(2)(i, j, k) = (z > z_ca && std::abs(y - y_ca) < 1.0 / 8) ? 0.1 : run.advectee(2)(i, j, k);

        run.advectee(3)(i, j, k) = 1 - 3. / 10 * (run.advectee(0)(i, j, k) + run.advectee(1)(i, j, k) + run.advectee(3)(i, j, k));

        run.g_factor()(i, j, k) = a * rho * cos(y);
      }
    }
  }
  
  solution = run.advectee();

  std::cout << "Calculating: " << filename << std::endl;

  run.advance(12 * day);

  T Linf = max(abs(run.advectee() - solution)) / max(abs(solution));
  T L1 = sum(run.g_factor() * abs(run.advectee() - solution)) / sum(run.g_factor() * abs(solution));
  T L2 = sqrt(sum(run.g_factor() * pow2(run.advectee() - solution)) / sum(run.g_factor() * pow2(solution)));

  std::cout << "real time:\t" << day << std::endl;
  std::cout << "Linf:\t" << Linf << std::endl;
  std::cout << "L1:\t" << L1 << std::endl;
  std::cout << "L2:\t" << L2 << std::endl;
}

int main()
{
  {
    enum { opts = opts::nug | opts::iga};
    const int opts_iters = 2;
    test<opts, opts_iters>("nug_iga", 90);
  }
}
