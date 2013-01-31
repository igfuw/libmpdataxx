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

  typedef typename homo_solver::mem_t::real_t real_t;

  private:

  virtual void forcings(real_t dt) = 0;

  protected:

  typedef homo_solver parent;

  // for use within forcings
  typename parent::mem_t::arr_t state(int e)
  {
    return this->mem.psi[e][this->mem.n[e]];
  }

  real_t dt;

  public:

  struct params : parent::params { real_t dt; };

  // 1D
  inhomo_solver(typename parent::mem_t &mem, const rng_t &i, params p) :
    parent(mem, i, p), dt(p.dt)
  {}

  // 2D
  inhomo_solver(typename parent::mem_t &mem, const rng_t &i, const rng_t &j, real_t dt) :
    parent(mem, i, j), dt(dt)
  {}

  // 3D
  // TODO

  void solve(int nt)
  {
    for (int t = 0; t < nt; ++t)
    {
      if (!naive) forcings(dt / 2);
      else forcings(dt);

      parent::solve(1);

      if (!naive) forcings(dt / 2);
    }
  }
};

template <class homo_solver>
using inhomo_solver_naive = inhomo_solver<homo_solver, true>;
