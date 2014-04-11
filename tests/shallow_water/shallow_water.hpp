/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/solvers/mpdata_rhs_vip.hpp> 
#include <libmpdata++/formulae/nabla_formulae.hpp>

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
template <typename ct_params_t, class enableif = void>
class shallow_water 
{};

template <class ct_params_t>
class shallow_water_common : public libmpdataxx::solvers::mpdata_rhs_vip<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata_rhs_vip<ct_params_t>;

  protected:

  // member fields
  const typename ct_params_t::real_t g;

  // 
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at);
    enum { n = 0 };    // just to make n, n+1 look nice :)
    assert(
      this->timestep == 0 && at == n 
      ||
      this->timestep  > 0 && at == n+1
    ); // note: we know only how to calculate R^{n+1}
       //       thus allowing to treat R^{n+1} as R^{n}
       //       only in the first timestep
  }

  void hook_post_step()
  {
    parent_t::hook_post_step();
    assert(min(this->psi_n(ct_params_t::ix::h)(this->ijk)) >= 0);  
  }

  void hook_ante_step()
  {
    parent_t::hook_ante_step();
    assert(min(this->psi_n(ct_params_t::ix::h)(this->ijk)) >= 0);  
  }

  public:

  // run-time parameters
  struct rt_params_t : parent_t::rt_params_t 
  {   
    typename parent_t::real_t g = 9.81; // default value 
  };

  // ctor
  shallow_water_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
    parent_t(args, p), 
    g(p.g)
  {}
};

// 1D version
template <typename ct_params_t>
class shallow_water<
  ct_params_t, 
  typename std::enable_if<ct_params_t::n_dims == 1>::type
> : public shallow_water_common<ct_params_t>
{
  static_assert(ct_params_t::n_eqns == 2, "{qx, h} in 1D");
  using parent_t = shallow_water_common<ct_params_t>;
  using parent_t::parent_t; // inheriting ctors
  using ix = typename ct_params_t::ix;

//<listing-1>
  void update_rhs(
    libmpdataxx::arrvec_t<
      typename parent_t::arr_t
    > &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    using namespace libmpdataxx::formulae::nabla;

    parent_t::update_rhs(rhs, dt, at);

    // due to grad() below -- TODO: include it into grad
    this->xchng(ix::h); 

    rhs.at(ix::qx)(this->i) -= 
      this->g 
      * this->psi_n(ix::h)(this->i) 
      * grad(this->psi_n(ix::h), this->i, this->di); 
  }
//</listing-1>
};

// 2D version
template <typename ct_params_t>
class shallow_water<
  ct_params_t, 
  typename std::enable_if<ct_params_t::n_dims == 2>::type
> : public shallow_water_common<ct_params_t>
{
  static_assert(ct_params_t::n_eqns == 3, "{qx, qy, h} in 2D");
  using parent_t = shallow_water_common<ct_params_t>;
  using parent_t::parent_t; // inheriting ctors
  using ix = typename ct_params_t::ix;

  template <int d, class arr_t>
  void forcings_helper(
    arr_t rhs, // TODO: ref?
    const libmpdataxx::rng_t &i,
    const libmpdataxx::rng_t &j,
    const typename ct_params_t::real_t &di
  )
  {
    using namespace libmpdataxx::formulae::nabla;
    rhs(pi<d>(i,j)) -= this->g * this->psi_n(ix::h)(pi<d>(i,j)) * grad<d>(this->psi_n(ix::h), i, j, di); 
  }

  /// @brief Shallow Water Equations: Momentum forcings for the X and Y coordinates
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    //
    parent_t::update_rhs(rhs, dt, at);

    //
    this->xchng(ix::h);

    //
    forcings_helper<0>(rhs.at(ix::qx), this->i, this->j, this->di);
    forcings_helper<1>(rhs.at(ix::qy), this->j, this->i, this->dj); 
  }
};
