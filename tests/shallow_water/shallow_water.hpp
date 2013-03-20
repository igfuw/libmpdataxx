/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <advoocat/solvers/detail/solver_velocity_common.hpp> // TODO: detail!
#include <advoocat/solvers/solver_inhomo.hpp> 

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
  * - \f$ \nabla_z = \partial_x + \partial_y \f$
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
  * momentum eq. plus u times mass continuity equation:
  * \f$ \partial_t (uh) + \nabla_z (uh) = -g h \nabla_z \eta \f$
  */
template <class solver_t, int qx, int qy, int h>
class shallow_water : public solvers::detail::solver_velocity_common<solvers::inhomo_solver<solver_t, solvers::strang>, qx, qy>
// TODO: check if qx, qy do not need some modification before being treated as u, w??
{
  using parent_t = solvers::detail::solver_velocity_common<solvers::inhomo_solver<solver_t, solvers::strang>, qx, qy>;

  public: 

  using real_t = typename parent_t::real_t;

  private:

  // member fields
  typename parent_t::arr_t dHdx, dHdy; // Blitz ctors have reference semantics - these will not be copies

  /// @brief Shallow Water Equations: Momentum forcings for the X and Y coordinates
  void forcings(real_t dt)  
  {
    // explicit formulation
    rng_t &i = this->i, &j = this->j;
 
    this->state(qx)(i,j) -= 0 + dHdx(i,j);
  }

  public:

  struct params_t : parent_t::params_t 
  {
    typename parent_t::arr_t dHdx, dHdy;
  };

  // ctor
  shallow_water( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) :
    parent_t(args, p), dHdx(p.dHdx), dHdy(p.dHdy)
  {}
};
