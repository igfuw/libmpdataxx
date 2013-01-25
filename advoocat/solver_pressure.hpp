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
//gradient, diveregnce
#include "nabla_formulae.hpp"

template <class inhomo_solver, int u, int w, int tht, int x, int z>
class pressure_solver : public inhomo_solver
{
  public:

  typedef typename inhomo_solver::real_t real_t;
  real_t ny;
  real_t nx;
  real_t dt;
  real_t Prs_amb;

  blitz::Array<real_t, 2> Prs;

  virtual void forcings(real_t dt) = 0;

  private:

  rng_t im, jm;
  //TODO don't assume dx=dz=1
  real_t dx = 1;
  real_t dz = 1;

  void extrp_velocity(int e) // extrapolate in time to t+1/2
  {            // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
    auto tmp = this->psi[e][this->n[e] -1];
    tmp /= -2;
    tmp += 3./2 * this->psi[e][this->n[e]]; 
  }

  void ini_courant()
  {
    this->xchng(u);      
    this->xchng(w);     

    courant::intrp<x>(this->C[x], this->psi[u][this->n[u]], im, this->j, this->dt, dx);
    courant::intrp<z>(this->C[z], this->psi[w][this->n[w]], jm, this->i, this->dt, dz);
  }

  void update_courant()
  {
    extrp_velocity(u);   //extrapolate velocity field in time (t+1/2)
    extrp_velocity(w);

    this->xchng(u, 1);      // filling halos for velocity filed
    this->xchng(w, 1);      // psi[n-1] was overwriten for that by extrp_velocity

    courant::intrp<x>(this->C[x], this->psi[u][this->n[u]-1], im, this->j, this->dt, dx);
    courant::intrp<z>(this->C[z], this->psi[w][this->n[w]-1], jm, this->i, this->dt, dz);
  } 

  void ini_pressure()
  {
    Prs(this->i, this->j) = this->Prs_amb;
  }

  void pressure_solver_update()
  {
    this->xchng(u);      
    this->xchng(w);     

    real_t beta = 1; //TODO

    this->bcx.fill_halos(this->Prs, this->j^this->halo);
    this->bcy.fill_halos(this->Prs, this->i^this->halo);

    this->Prs += beta * nabla_op::div(
      //TODO add density field
      /*real_t(1) */ this->psi[u] - dt/2 * nabla_op::grad<0>(Prs, this->i, this->j, real_t(1)),
      /*real_t(1) */ this->psi[w] - dt/2 * nabla_op::grad<1>(Prs, this->j, this->i, real_t(1)),
      this->i,
      this->j,
      real_t(1),   //dx  TODO
      real_t(1)    //dy
    );

//    tmp_u = - nabla_op::grad<0>(p-prs_amb/rho, this->i, this->j, real_t(1));
//    tmp_w = - nabla_op::grad<1>(p-prs_amb/rho, this->j, this->i, real_t(1));

  }


  public:
  // ctor
  pressure_solver(int nx, int ny, real_t dt, real_t Prs_amb) :
    inhomo_solver(nx, ny, dt), dt(dt),
    ny(ny), nx(nx),
    Prs_amb(Prs_amb),
    im(this->i.first() - 1, this->i.last()),
    jm(this->j.first() - 1, this->j.last())
  {}

  void solve(int nt)
  {
    for (int t = 0; t < nt; ++t)
    {
      if (t==0)
      {
        ini_courant();
        ini_pressure();
        forcings(dt / 2);
        inhomo_solver::parent::solve(1);
        forcings(dt / 2);
        pressure_solver_update();
//pressure_solver_apply
      }
      if (t!=0)
      {
        update_courant();
        forcings(dt / 2);
//pressure_solver_apply
        inhomo_solver::parent::solve(1);
        // this->xchng_all();
        // this->advop_all();
        // this->cycle_all(); 
        forcings(dt / 2);
        pressure_solver_update();
//pressure_solver_apply
      }
    }
  }
};
