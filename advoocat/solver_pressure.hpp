/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include "solver_common.hpp"

template <class homo_solver>
class pressure_solver : public homo_solver
{
  //TODO as a template
  enum{tht, prs, u, w};

  using real_t = typename homo_solver::real_t;

  //arrvec_t<typename homo_solver::arr_t> rhs;
  real_t dt;
  virtual void forcings(real_t dt) = 0;

  void update_courant(real_t dt)
  {
    // psi[n-1] will not be used anymore, and
    // it will be intentionally overwritten!
    arr_2d_t<real_t> tmp;
    tmp = this->psi[u][n-1];
    tmp /= -2;
    tmp += 3./2 * this->psi[u][n];
    //this->C = interpol(tmp); 

    const rng_t &i = this->i;
    const rng_t &j = this->j;

    //filing halos for velocity filed
    //TODO remeber to take those that were extrapolated before
    this->xchng(u);
    this->xchng(w);
  
    //TODO - don't assume dx=1
    this->courant(0)(i+h,j)   = .5 / dt * (this->state(u)(i^h,j  )+this->state(u)((i^h)+1, j     ));
    this->courant(1)(i  ,j+h) = .5 / dt * (this->state(w)(i  ,j^h)+this->state(w)( i     ,(j^h)+1));
  }

  public:
  // ctor
  pressure_solver(int nx, int ny, real_t dt) :
    homo_solver(nx, ny), dt(dt)
  {}

  void solve(int nt)
  {
    for (int t = 0; t < nt; ++t)
    {
      //extrapolate velocity field in time (t+1/2)
      //interpolate velocity to fill courant field
      update_courant(dt);

      forcings(dt / 2);
      this->xchng_all();
      this->advop(tht);
      this->advop(u); 
      this->advop(w);   
      this->cycle_all(); 
      forcings(dt / 2);
    }
  }
};

