/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

/// @brief tha main namescape
namespace libmpdataxx
{};

/** @mainpage
  *
  * @section INTRODUCTION
  *
  * libmpdata++ is C++ library of parallel, object-oriented implementations
  * of the MPDATA family solvers of generalised transport equations 
  * of the form:
  *
  * \f$ \partial_t \psi + \nabla \cdot (\vec{v} \psi) = R \f$
  *
  * where
  *
  * \f$ \psi = [\psi_1, \psi_2, \ldots ] \f$ is a set of conservative dependent variables, 
  * \f$ R = [R_1, R_2, \ldots ] \f$ are the forcing terms, 
  * and \f$ \vec{v} = [u, v, w] \f$ is the velocity field.
  *
  * The theory behind MPDATA solvers was developed by Piotr Smolarkiewicz et al.
  * (see e.g. @copydetails bib::Smolarkiewicz_2006, for a review and list of references).
  * Development of libmpdata++ is carried out by Sylwester Arabas, Anna Jaruga and
  * co-workers at the [Institute of Geophysics](http://www.igf.fuw.edu.pl/),
  * [Faculty of Physics](http://www.fuw.edu.pl/),
  * [University of Warsaw](http://www.uw.edu.pl/) (the copyright holder)
  * with funding from the Polish [National Science Centre](http://www.ncn.gov.pl/).
  *
  * libmpdata++ is based on the Blitz++ library for high-performance array handling.
  * libmpdata++ (and Blitz++) is a header-only library.
  * 
  * @section sec_concepts KEY CONCEPTS
  *
  * Shortest example:
  * TODO
  * 
  * Library "layers":
  * - solvers
  * - boundary conditions
  * - concurency schemes
  * - output mechanisms
  *
  * At present libmpdata++ handles 1D, 2D and 3D computations in cartesian coordinates on an Arakawa-C grid.
  *
  * @section sec_examples SUMMARY OF SELECTED EXAMPLES
  *
  * | test name                                               | n_dims | n_eqs | velocity       | output                                      | concurrency                                  | system solved                                        | sample figure (svg)                                                             |
  * | :-----------------------------------------------------: | :----: | :---: | :------------: | :-----------------------------------------: | :------------------------------------------: | :--------------------------------------------------: | :-----------------------------------------------------------------------------: |
  * | \ref test_gnuplot-iostream_1d.cpp "gnuplot-iostream_1d" | 1      | 1     | constant       | \ref libmpdataxx::output::gnuplot "gnuplot" | \ref libmpdataxx::concurr::threads "threads" | \f$ \partial_t \psi + \nabla (\vec{u} \psi) = 0 \f$  | \image html "../../tests/gnuplot-iostream_1d/figure_iters=2.svg"                |
  * | \ref test_gnuplot-iostream_2d.cpp "gnuplot-iostream_2d" | 2      | 1     | constant       | \ref libmpdataxx::output::gnuplot "gnuplot" | \ref threads.hpp "threads (openmp/boost)"    | \f$ \partial_t \psi + \nabla (\vec{u} \psi) = 0 \f$  | \image html "../../tests/gnuplot-iostream_2d/figure_iters=2_psi_96.svg"         |
  * | \ref test_rotating_cone.cpp "rotating_cone"             | 2      | 1     | prescribed     | \ref libmpdataxx::output::gnuplot "gnuplot" | \ref threads.hpp "threads (openmp/boost)"    | \f$ \partial_t \psi + \nabla (\vec{u} \psi) = 0 \f$  | \image html "../../tests/rotating_cone/figure_iters=2_psi_0.svg"                |
  * | \ref test_harmosc.cpp "harmosc"                         | 1      | 2     | constant       | \ref libmpdataxx::output::gnuplot "gnuplot" | \ref threads.hpp "threads (openmp/boost)"    | pair of coupled 1D harmonic oscillators: \f$ \\ \partial_t \psi + \nabla (\vec{u} \psi) =  \omega \phi \\ \partial_t \phi + \nabla (\vec{u} \phi) = -\omega \psi \f$ | \image html "../../tests/harmosc/figure_euler_it=1.svg" |
  * | \ref test_shallow_water.cpp "shallow_water"             | 2      | 3     | resolved       | \ref libmpdataxx::output::gnuplot "gnuplot" | \ref threads.hpp "threads (openmp/boost)"    | 2D shallow-water system: \f$ \\ \partial_t h + \nabla_z (h \cdot u) = 0 \\ \partial_t (uh) + \nabla (\vec u \cdot uh) = -g h \partial_x h \\ \partial_t (vh) + \nabla (\vec u \cdot vh) = -g h \partial_y h \f$ | \image html "../../tests/shallow_water/figure_h_100.svg" |
  * | \ref test_isentropic.cpp "isentropic"                   | 1      |       | resolved       | \ref libmpdataxx::output::netcdf "netcdf"   | \ref threads.hpp "threads (openmp/boost)"    | TODO | TODO |
  * | \ref test_bombel.cpp "bombel"                           | 2      | 3     | resolved       | \ref libmpdataxx::output::gnuplot "gnuplot" | \ref threads.hpp "threads (openmp/boost)"    | 2D Navier-Stokes (adiabatic, constant density): \f$ \\ \partial_t \theta + \nabla (\vec{u} \theta) = 0 \\ \partial_t u + \nabla (\vec{u} u) = -\partial_x \frac{p-p_e}{\rho} \\ \partial_t w + \nabla (\vec{u} w) = -\partial_z \frac{p-p_e}{\rho} + g \frac{\theta - \theta_e}{\theta_e} \f$ | \image html "../../tests/bombel/figure_pc_u_10.svg"|
  *
  *
  *
  * @section sec_install INSTALLATION AND USAGE
  * suggested compiler options (by compiler): -march=native, -Ofast, -DNDEBUG, -lblitz (opt), -DBZDEBUG (opt), -std=c++11
  */
// TODO: mention:
 /*
  *   - Richardson solver 
  *   - Minimum-Residual solver (solver_pressure_mr.hpp)
  *   - Conjugate-Residual solver (solver_pressure_cr.hpp)
  *   - Conjugate-Residual with preconditioner (solver_pressure_pc.hpp)
  */

/** @page HACKING HACKING file (coding conventions)
 *  @section sec_HACKING HACKING file (coding conventions)
 *  @verbinclude "../HACKING"
 */
/** @page COPYING GNU General Public License version 3
 *  @section sec_COPYING GNU General Public License version 3
 *  @verbinclude "../COPYING"
 */
