/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief a central reference list for the documentation
  */

#error "intended for use with Doxygen, not to be included during compilation"

/// @brief bibliography items (trick to emulate BiBTeX-like logic in Doxygen)
namespace bib
{
  /// @brief
  /// [Anderson and Fattahi 1974]: http://dx.doi.org/10.1175/1520-0469(1974)031<1500:ACONSO>2.0.CO;2
  /// [Anderson and Fattahi 1974][]
  /// @details
  /// A comparison of numerical solutions of the advective equation.
  /// J. Atmos. Sci., 31, 1500-1506.
  typedef void Anderson_and_Fattahi_1974;

  /// @brief
  /// [Rasiński et al. 2011]: http://dx.doi.org/10.1016/j.atmosres.2011.06.020
  /// [Rasiński et al. 2011][]
  /// @details
  /// Observations and kinematic modeling of drizzling marine stratocumulus.
  /// Atmos. Res., 102, 120-135.
  typedef void Rasinski_et_al_2011;

  /// @brief
  /// [Smolarkiewicz 1984]: http://dx.doi.org/10.1016/0021-9991(84)90121-9
  /// [Smolarkiewicz 1984][]
  /// @details
  /// A fully multidimensional positive definite advection transport algorithm with small implicit diffusion
  /// J. Comp. Phys., 54, 352-362.
  typedef void Smolarkiewicz_1984;

  /// @brief
  /// [Smolarkiewicz 2006]: http://dx.doi.org/10.1002/fld.1071
  /// [Smolarkiewicz 2006][]
  /// @details
  /// Multidimensional positive deÿnite advection transport algorithm: An overview
  /// Intl. J. Numer. Meth. Fluids., 50, 1123-1144.
  typedef void Smolarkiewicz_2006;

  /// @brief 
  /// [Smolarkiewicz_and_Grabowski_1990]: http://dx.doi.org/TODO
  /// [Smolarkiewicz_and_Grabowski_1990][]
  /// @details
  /// The Multidimensional Positive Definite Advection Transport Algorithm: Nonoscillatory Option
  /// J. Comp. Phys., 86, 355-375.
  typedef void Smolarkiewicz_and_Grabowski_1990;  

}; // namespace bib
