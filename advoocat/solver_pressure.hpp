/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include "solver_common.hpp"
#include "solver_inhomo.hpp"
#include "courant_formulae.hpp"

template <class inhomo_solver, int u, int w, int tht, int x, int z>
class pressure_solver : public inhomo_solver
{
  public:

  typedef typename inhomo_solver::real_t real_t;
  real_t ny;
  real_t nx;

  private:

  rng_t im, jm;
  //TODO don't assume dx=dz=1
  real_t dx = 1;
  real_t dz = 1;

  void extrp_velocity(int e, int t) //extrapolate in time to t+1/2
  {                  
    // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
    auto tmp = this->psi[e][this->n[e] -1];
    if(t == 0)
    {
      tmp = this->psi[e][this->n[e]];      
    }
    if(t != 0)
    {
      tmp /= -2;
      tmp += 3./2 * this->psi[e][this->n[e]]; 
    }
  }

  void update_courant(int t)
  {
    extrp_velocity(u, t);   //extrapolate velocity field in time (t+1/2)
    extrp_velocity(w, t);

    this->xchng(u, 1);      // filling halos for velocity filed
    this->xchng(w, 1);      // psi[n-1] was overwriten for that by extrp_velocity

    courant::intrp<x>(this->C[x], this->psi[u][this->n[u]-1], im, this->j, this->dt, dx);
    courant::intrp<z>(this->C[z], this->psi[w][this->n[w]-1], jm, this->i, this->dt, dz);
  } 

  public:
  // ctor
  pressure_solver(int nx, int ny, real_t dt) :
    inhomo_solver(nx, ny, dt),
    ny(ny), nx(nx),
    im(this->i.first() - 1, this->i.last()),
    jm(this->j.first() - 1, this->j.last())
  {}

  void solve(int nt)
  {
    for (int t = 0; t < nt; ++t)
    {
      update_courant(t);
      inhomo_solver::solve(1);
      // forcings(dt / 2);
      // this->xchng_all();
      // this->advop_all();
      // this->cycle_all(); 
      // forcings(dt / 2);
    }
  }
};
