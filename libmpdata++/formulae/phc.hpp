/** @file
 *  author Anna Jaruga
 *  @copyright University of Warsaw
 *  @date January 2013
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the phc class encapsulating 
 *    a catalogue of physical constants
 */
#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/units.hpp>

// TODO: remove this file completely

#define phc_decltype_return(expr) -> decltype(expr) { return expr; }

#define phc_declare_const_macro(name, value, unit) template <typename real_t> \
  static constexpr auto name() phc_decltype_return(real_t(value) * unit)

#define phc_derived_const_macro(name, value) template <typename real_t> \
  static constexpr auto name() phc_decltype_return(value)

#define phc_declare_funct_macro template <typename real_t> 

#include <boost/math/constants/constants.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    // pi
    phc_derived_const_macro(pi, boost::math::constants::pi<real_t>())

    // acceleration due to gravity
    phc_declare_const_macro(g, 9.81, si::metres_per_second_squared)

    // pressure in the definition of potential temperature
    phc_declare_const_macro(p_1000, 100000, si::pascals)

    // molar masses
    phc_declare_const_macro(M_d, 0.02896, si::kilograms / si::moles) // dry air

    // universal gas constant (i.e. the Boltzmann times the Avogadro constants)
    phc_declare_const_macro(kaBoNA, 8.314472, si::joules / si::kelvins / si::moles)

    // gas constant for dry air
    phc_derived_const_macro(R_d, kaBoNA<real_t>() / M_d<real_t>())

    // specific heat capacities
    phc_declare_const_macro(c_pd, 1005, si::joules / si::kilograms / si::kelvins) // dry air

    // Exner function exponent for dry air
    phc_derived_const_macro(R_d_over_c_pd, R_d<real_t>() / c_pd<real_t>())

  }; // namespace formulae
}; // namespace libmpdataxx
