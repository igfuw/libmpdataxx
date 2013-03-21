/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <advoocat/solvers/detail/solver_velocity_common.hpp> // TODO: bez detail!
#include <advoocat/solvers/solver_inhomo.hpp> 
#include <advoocat/formulae/nabla_formulae.hpp>
#include <advoocat/formulae/phc.hpp>

/** @brief the 2D shallow-water equations system
  *
  * Consult chapter 3 in Vallis 2008 for a detailed derivation.
  *
  * The key assumptions are:
  * - horizontal scale is much larger than the vertical scale (\f$ u \approx u(x) \f$)
  * - hydrostatic equillibrium
  * - constant density
  *
  * Nomenclature:
  * - \f$ \eta(x,y) \f$ - (absolute) height of the fluid surface
  * - \f$ \eta_0(x,y) \f$ - bathymetry
  * - \f$ h = \eta - \eta_0 \f$ - thickness of the fluid layer
  * - \f$ \vec{u} = (u,v) \f$
  * - \f$ \nabla_z = (\partial_x, \partial_y) \f$ 
  *
  * momentum equation:
  * \f$ \partial_t u + u \cdot \nabla_z u = - \frac{1}{\rho} \nabla_z p \f$
  *
  * pressure in a column of the constant-density fluid:
  * \f$ p = p_0 - \rho g z = p_0 + \rho g \cdot (\eta(x) - z) \f$
  *
  * mass continuity equation:
  * \f$ \partial_t h + \nabla_z (h \cdot u) = 0 \f$
  *
  * h times momentum eq. plus u times mass continuity equation:
  * \f$ \partial_t (uh) + \nabla_z (u \cdot uh) = -g h \nabla_z \eta \f$
  */
template <typename real_t, int n_iters, int qx, int qy, int h>
class shallow_water : public solvers::detail::solver_velocity_common<
  solvers::inhomo_solver<
    solvers::mpdata_2d<real_t, n_iters, 3>, 
    solvers::strang
  >, qx, qy, h
>
{
  using parent_t = solvers::detail::solver_velocity_common<
    solvers::inhomo_solver<
      solvers::mpdata_2d<real_t, n_iters, 3>, 
      solvers::strang
    >, qx, qy, h
  >;

  private:

  const real_t g;

  template <int d, class arr_t>
  void forcings_helper(
    arr_t qq,
    const arr_t hh,
    const rng_t &i,
    const rng_t &j,
    const real_t &dx,
    const real_t &dt
  )
  {
    using namespace formulae::nabla;
    qq(pi<d>(i,j)) -= dt * g * hh(pi<d>(i,j)) * grad<d>(hh, i, j, dx); 
  }

  /// @brief Shallow Water Equations: Momentum forcings for the X and Y coordinates
  void forcings(real_t dt)  
  {
    this->xchng(h);
    forcings_helper<0>(this->state(qx), this->state(h), this->i, this->j, this->dx, dt);
    forcings_helper<1>(this->state(qy), this->state(h), this->j, this->i, this->dz, dt); // TODO: rename dz->dy?
  }

  public:

  // ctor
  shallow_water( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::params_t &p
  ) :
    parent_t(args, p), 
    g(formulae::g<real_t>() * si::seconds * si::seconds / si::metres)
  {}
};
