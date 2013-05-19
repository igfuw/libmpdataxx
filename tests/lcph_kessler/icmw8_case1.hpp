#pragma once

// TODO: includes!

// 8th ICMW case 1 by Wojciech Grabowski)
namespace icmw8_case1
{
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta = libcloudphxx::common::theta_std;

  const quantity<si::temperature, real_t> 
    th_0 = 289 * si::kelvins;
  const quantity<si::dimensionless, real_t> 
    rv_0 = 7.5e-3;
  const quantity<si::pressure, real_t> 
    p_0 = 101500 * si::pascals;
  const quantity<si::velocity, real_t> 
    w_max = real_t(.6) * si::metres_per_second;
  const quantity<si::length, real_t> 
    z_0  = 0    * si::metres,
    nzdz = 1500 * si::metres, 
    nxdx = 1500 * si::metres;
  const quantity<si::time, real_t>
    dt = 4 * si::seconds;

  // density profile as a function of altitude
  real_t rhod(real_t z)
  {
    quantity<si::pressure, real_t> p = hydrostatic::p(
      z * si::metres, th_0, rv_0, z_0, p_0
    );
    
    quantity<si::mass_density, real_t> rhod = theta::rhod(
      p, th_0, rv_0
    );

    return rhod / si::kilograms * si::cubic_metres;
  }

  // to make the rhod() function accept Blitz arrays as arguments
  BZ_DECLARE_FUNCTION(rhod);


  /// (similar to eq. 2 in @copydetails Rasinski_et_al_2011, Atmos. Res. 102)
  /// @arg xX = x / (nx*dx)
  /// @arg zZ = z / (nz*dz)
  real_t psi(real_t xX, real_t zZ)
  {
    return - sin(pi<real_t>() * zZ) * cos(2 * pi<real_t>() * xX);
  }
  BZ_DECLARE_FUNCTION2_RET(psi, real_t)

  // function expecting a libmpdata solver parameters struct as argument
  template <class T>
  void setopts(T &params, int nz)
  {
    params.dt = dt / si::seconds;
    params.dz = nzdz / si::metres / nz;

    params.rhod.resize(nz);
    for (int j = 0; j < nz; ++j)
      params.rhod[j] = rhod((j+.5) * params.dz);
  }


  // function expecting a libmpdata++ solver as argument
  template <class T>
  void intcond(T &solver)
  {
    // helper ondex placeholders
    blitz::firstIndex i;
    blitz::secondIndex j;

    // dx, dy ensuring 1500x1500 domain
    int 
      nx = solver.state().extent(x),
      nz = solver.state().extent(z);
    real_t 
      dx = nxdx / si::metres / nx, 
      dz = nzdz / si::metres / nz; 

    // constant potential temperature & water vapour mixing ratio profiles
    solver.state(rhod_th_ix) = rhod((j+.5)*dz) * (th_0 / si::kelvins); // TODO: should be theta_dry and is theta
    solver.state(rhod_rv_ix) = rhod((j+.5)*dz) * real_t(rv_0);

    // TODO: only for bulk!
    solver.state(rhod_rc_ix) = 0;
    solver.state(rhod_rr_ix) = 0;

    // velocity field obtained by numerically differentiating a stream function
    {
      real_t A = (w_max / si::metres_per_second) * solver.state().extent(x) * dx / pi<real_t>();
    
      solver.courant(x) = - A * (
	psi(i/real_t(nx), (j+.5+.5)/nz)-
	psi(i/real_t(nx), (j+.5-.5)/nz)
      ) / dz             // numerical derivative
      / rhod((j+.5)* dz) // psi defines rho_d times velocity
      * (dt / si::seconds) / dx;         // converting to Courant number

      solver.courant(z) = A * (
	psi((i+.5+.5)/nx, j/real_t(nz)) -
	psi((i+.5-.5)/nx, j/real_t(nz))
      ) / dx 
      / rhod(j * dz)
      * (dt / si::seconds) / dz; 
    }
  }
};
