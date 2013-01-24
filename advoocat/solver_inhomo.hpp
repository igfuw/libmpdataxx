/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief improved Euler inhomogeneous solver  
 *        (cf. eq. 32 in Smolarkiewicz 1998)
 */

#pragma once

template <class homo_solver, bool naive = false>
class inhomo_solver : public homo_solver
{
  public:

  typedef typename homo_solver::real_t real_t;

  private:

  virtual void forcings(real_t dt) = 0;

  protected:

  real_t dt;

  public:

  // 1D
  inhomo_solver(int n, real_t dt) :
    homo_solver(n), dt(dt)
  {}

  // 2D
  inhomo_solver(int nx, int ny, real_t dt) :
    homo_solver(nx, ny), dt(dt)
  {}

  void solve(int nt)
  {
    for (int t = 0; t < nt; ++t)
    {
      if (!naive) forcings(dt / 2);
      else forcings(dt);

      homo_solver::solve(1);

      if (!naive) forcings(dt / 2);
    }
  }
};

template <class homo_solver>
using inhomo_solver_naive = inhomo_solver<homo_solver, true>;
