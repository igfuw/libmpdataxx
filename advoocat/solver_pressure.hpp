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

  blitz::Array<real_t, 2> Phi;
  //TODO probably don't need those
  blitz::Array<real_t, 2> tmp_u;
  blitz::Array<real_t, 2> tmp_w;
  blitz::Array<real_t, 2> tmp_div;
  blitz::Array<real_t, 2> tmp_x;
  blitz::Array<real_t, 2> tmp_z;

  virtual void forcings(real_t dt) = 0;

  private:

  rng_t im, jm, i_domain, j_domain;
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

    courant::intrp<x>(this->C[x], this->psi[u][this->n[u]], im, this->j^this->halo, this->dt, dx);
    courant::intrp<z>(this->C[z], this->psi[w][this->n[w]], jm, this->i^this->halo, this->dt, dz);
  }

  void update_courant()
  {
    extrp_velocity(u);      //extrapolate velocity field in time (t+1/2)
    extrp_velocity(w);

    this->xchng(u, 1);      // filling halos for velocity filed
    this->xchng(w, 1);      // psi[n-1] was overwriten for that by extrp_velocity

    courant::intrp<x>(this->C[x], this->psi[u][this->n[u]-1], im, this->j^this->halo, this->dt, dx);
    courant::intrp<z>(this->C[z], this->psi[w][this->n[w]-1], jm, this->i^this->halo, this->dt, dz);
  } 

  void ini_pressure()
  {
    // dt/2 * (Prs-Prs_amb) / rho
    Phi(this->i^this->halo, this->j^this->halo) = real_t(0);
  }

  void pressure_solver_update(real_t dt)
  {
    real_t beta = .25;  //TODO
    real_t rho = 1.;   //TODO    

    int halo = this->halo;
    rng_t i = this->i;
    rng_t j = this->j;

    tmp_u = this->psi[u][this->n[u]];  //making copy of velocity field before entering pressure solver 
    tmp_w = this->psi[w][this->n[w]];  //TODO not needed?

std::cerr<<"--------------------------------------------------------------"<<std::endl;
    //pseudo-time loop
    real_t err = 1.;
    while (err > .01)
    {
      this->xchng(Phi,   i^halo, j^halo);
      this->xchng(tmp_u, i^halo, j^halo);
      this->xchng(tmp_w, i^halo, j^halo);

      tmp_x(i, j) = rho * tmp_u(i, j) - nabla_op::grad<0>(Phi(i^halo, j^halo), i, j, real_t(1));
      tmp_z(i, j) = rho * tmp_w(i, j) - nabla_op::grad<1>(Phi(i^halo, j^halo), j, i, real_t(1));
 
      this->xchng(tmp_x, i^halo, j^halo);
      this->xchng(tmp_z, i^halo, j^halo);

      tmp_div(i, j) = nabla_op::div(tmp_x(i^halo,j^halo), tmp_z(i^halo, j^halo), i, j, real_t(1), real_t(1));
      Phi(i, j) -= beta * tmp_div(i, j);

      this->xchng(tmp_u, i^halo, j^halo);
      this->xchng(tmp_w, i^halo, j^halo);

      err = std::abs(max(tmp_div)) + std::abs(min(tmp_div));
std::cerr<<"div:  ( "<<min(tmp_div)<<" --> "<<max(tmp_div)<<" )"<<std::endl;
    }
    //end of pseudo_time loop
    this->xchng(this->Phi, i^halo, j^halo);

    tmp_u(i, j) -= nabla_op::grad<0>(Phi(i^halo, j^halo), i, j, real_t(1));
    tmp_w(i, j) -= nabla_op::grad<1>(Phi(i^halo, j^halo), j, i, real_t(1));

    tmp_u -= this->psi[u][this->n[u]];
    tmp_w -= this->psi[w][this->n[w]];
  }

  void pressure_solver_apply(real_t dt)
  {
    auto U = this->psi[u][this->n[u]];
    auto W = this->psi[w][this->n[w]];

    U += tmp_u;
    W += tmp_w;
  }

  public:
  // ctor
  pressure_solver(int nx, int ny, real_t dt, real_t Prs_amb) :
    inhomo_solver(nx, ny, dt), dt(dt),
    ny(ny), nx(nx),
    Prs_amb(Prs_amb),
    Phi(this->i^this->halo, this->j^this->halo), //perturbation of pressure for pressure solver
    tmp_x(this->i^this->halo, this->j^this->halo),
    tmp_z(this->i^this->halo, this->j^this->halo),
    tmp_div(this->i, this->j),
    tmp_u(this->i^this->halo, this->j^this->halo),
    tmp_w(this->i^this->halo, this->j^this->halo),
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
        pressure_solver_update(dt);
        pressure_solver_apply(dt);
      }
      if (t!=0)
      {
std::cerr<<"t= "<<t<<std::endl;
        update_courant();
        forcings(dt / 2);
        pressure_solver_apply(dt);
        inhomo_solver::parent::solve(1);
        // this->xchng_all();
        // this->advop_all();
        // this->cycle_all(); 
        forcings(dt / 2);
        pressure_solver_update(dt);
        pressure_solver_apply(dt);
      }
    }
  }
};
