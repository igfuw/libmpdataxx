/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief a central reference list for the documentation
  */

/// @brief bibliography items (trick to emulate BiBTeX-like logic in Doxygen)
namespace bib
{
  /// @brief
  /// A comparison of numerical solutions of the advective equation.
  /// J. Atmos. Sci., 31, 1500-1506.
  /// @details
  /// [Anderson and Fattahi 1974]: http://dx.doi.org/10.1175/1520-0469(1974)031<1500:ACONSO>2.0.CO;2
  /// [Anderson and Fattahi 1974][]
  typedef void Anderson_and_Fattahi_1974;

  /// @brief
  /// A fully multidimensional positive definite advection transport algorithm with small implicit diffusion
  /// J. Comp. Phys., 54, 352-362.
  /// @details
  /// [Smolarkiewicz 1984]: http://dx.doi.org/10.1016/0021-9991(84)90121-9
  /// [Smolarkiewicz 1984][]
  typedef void Smolarkiewicz_1984;

  /// @brief
  /// Multidimensional positive deÿnite advection transport algorithm: An overview
  /// Intl. J. Numer. Meth. Fluids., 50, 1123-1144.
  /// @details
  /// [Smolarkiewicz 2006]: http://dx.doi.org/10.1002/fld.1071
  /// [Smolarkiewicz 2006][]
  typedef void Smolarkiewicz_2006;
}; // namespace bib
