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
  //TODO probably don't need those
  blitz::Array<real_t, 2> tmp_u;
  blitz::Array<real_t, 2> tmp_w;
  blitz::Array<real_t, 2> tmp;
  blitz::Array<real_t, 2> tmptmp;
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
    Prs(this->i^this->halo, this->j^this->halo) = this->Prs_amb;
  }

  void pressure_solver_update(real_t dt)
  {
std::cerr<<"------------pressure solver update-------------"<<std::endl;
 
    real_t beta = .5; //TODO
    real_t rho = 1.;  //TODO    

    this->xchng(u);      
    this->xchng(w);     
    this->bcx.fill_halos(this->Prs, this->j^this->halo);
    this->bcy.fill_halos(this->Prs, this->i^this->halo);

std::cerr<<"pressure before solver"<<std::endl;
std::cerr<<"Prs = "<<Prs(this->i, this->j)<<std::endl;

    tmp_u = this->psi[u][this->n[u]];  //making copy of velocity field before entering pressure solver 
    tmp_w = this->psi[w][this->n[w]];  //TODO not needed?
//---------------pseudo-time loop
std::cerr<<"pseudo_time loop:"<<std::endl;
    for (int tau = 0; tau <= 10; tau++)
    {
std::cerr<<"tau= "<<tau<<std::endl;
tmp_x(this->i, this->j) = rho * (tmp_u - dt/2 * nabla_op::grad<0>(Prs(this->i^this->halo, this->j^this->halo), this->i, this->j, real_t(1)));
tmp_z(this->i, this->j) = rho * (tmp_w - dt/2 * nabla_op::grad<1>(Prs(this->i^this->halo, this->j^this->halo), this->j, this->i, real_t(1)));
 
this->bcx.fill_halos(tmp_x, this->j^this->halo);
this->bcy.fill_halos(tmp_x, this->i^this->halo);
this->bcx.fill_halos(tmp_z, this->j^this->halo);
this->bcy.fill_halos(tmp_z, this->i^this->halo);

//if (tau == 10) std::cerr<<"rho * (u - grad_x(p)) "<<tmp_x(this->i^this->halo, this->j^this->halo)<<std::endl;
//if (tau == 10) std::cerr<<"rho * (w - grad_z(p)) "<<tmp_z(this->i^this->halo, this->j^this->halo)<<std::endl;

tmptmp(this->i, this->j) = nabla_op::div(
                             tmp_x(this->i^this->halo, this->j^this->halo), 
                             tmp_z(this->i^this->halo, this->j^this->halo),
                             this->i,this->j, real_t(1), real_t(1));

if (tau ==10) std::cerr<<"div(rho * (u - grad(p))) "<<tmptmp<<std::endl;

      Prs(this->i, this->j) += beta * tmptmp(this->i, this->j);
//        nabla_op::div(
//        //TODO add density field
//        tmp_x(this->i^this->halo, this->j^this->halo),
//        tmp_z(this->i^this->halo, this->j^this->halo),
//        this->i,
//        this->j,
//        real_t(1),   //dx  TODO
//        real_t(1)    //dy
//      );

if (tau ==10) std::cerr<<"pressure after solver"<<std::endl;
if (tau ==10) std::cerr<<"(should be different than initial pressure)"<<std::endl;
if (tau ==10) std::cerr<<"Prs = Prs + div(rho * (u -grad(p)) "<<Prs(this->i, this->j)<<std::endl;

      this->bcx.fill_halos(this->Prs, this->j^this->halo);
      this->bcy.fill_halos(this->Prs, this->i^this->halo);

      tmp_u(this->i, this->j) -= nabla_op::grad<0>(Prs(this->i^this->halo, this->j^this->halo)-Prs_amb/rho, this->i, this->j, real_t(1));
      tmp_w(this->i, this->j) -= nabla_op::grad<1>(Prs(this->i^this->halo, this->j^this->halo)-Prs_amb/rho, this->j, this->i, real_t(1));
    }
std::cerr<<"end of pseudo_time loop"<<std::endl;

//----------------------------------------------
    tmp_u(this->i, this->j) -= this->psi[u][this->n[u]];
    tmp_w(this->i, this->j) -= this->psi[w][this->n[w]];
  }

  void pressure_solver_apply(real_t dt)
  {
    auto U = this->psi[u][this->n[u]];
    auto W = this->psi[w][this->n[w]];

    U += tmp_u * dt/2.;
    W += tmp_w * dt/2.;
  }

  public:
  // ctor
  pressure_solver(int nx, int ny, real_t dt, real_t Prs_amb) :
    inhomo_solver(nx, ny, dt), dt(dt),
    ny(ny), nx(nx),
    Prs_amb(Prs_amb),
    Prs(this->i^this->halo, this->j^this->halo),
    tmp_u(this->i^this->halo, this->j^this->halo),
    tmp_w(this->i^this->halo, this->j^this->halo),
    tmp_x(this->i^this->halo, this->j^this->halo),
    tmp_z(this->i^this->halo, this->j^this->halo),
    tmp(this->i, this->j),
    tmptmp(this->i, this->j),
    im(this->i.first() - 1, this->i.last()),
    jm(this->j.first() - 1, this->j.last())
  {}

  void solve(int nt)
  {
    for (int t = 0; t < nt; ++t)
    {
      if (t==0)
      {
std::cerr<<"initial step"<<std::endl;
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
std::cerr<<"timestepping t= "<<t<<std::endl;
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
