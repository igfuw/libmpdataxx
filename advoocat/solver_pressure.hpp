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
  // u, w must be first and second
  enum{u, w, tht, prs};

  using real_t = typename homo_solver::real_t;

  real_t dt;
  virtual void forcings(real_t dt) = 0;

  void extrp_velocity(int e, int t) //extrapolate in time to t+1/2
  {                  // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
    auto tmp = this->psi[e][this->n[e] + (t!=0?-1:0)];
    tmp /= -2;
    tmp += 3./2 * this->psi[e][this->n[e]]; 
  }

  void intrp_courant(real_t dt) //interpolate in space to courant field
  { 
    const rng_t &i = this->i;
    const rng_t &j = this->j;

    this->xchng(u, 1);  // filling halos for velocity filed
    this->xchng(w, 1);  // psi[n-1] was overwriten for that by extrp_velocity
                        // so in xchng previous time level has to be used

    auto tmp_u = this->psi[u][this->n[u]-1];
    auto tmp_w = this->psi[w][this->n[w]-1];
  
    //TODO - don't assume dx=1
    real_t dx = 1;
    this->courant(u)(i+h,j  ) = dt / dx * .5 * (tmp_u(i^h,j  )+tmp_u((i^h)+1, j     ));
    this->courant(w)(i  ,j+h) = dt / dx * .5 * (tmp_w(i  ,j^h)+tmp_w( i     ,(j^h)+1));
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
      extrp_velocity(u, t);   //extrapolate velocity field in time (t+1/2)
      extrp_velocity(w, t);
      intrp_courant(dt);      //interpolate velocity to fill courant field
 
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

