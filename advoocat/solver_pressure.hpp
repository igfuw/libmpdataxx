/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once
#include "solvers_common.hpp"

template <class homo_solver>
class pressure_solver : public homo_solver
{
  //TODO as a template
  enum{tht, prs, u, w};

  //arrvec_t<typename homo_solver::arr_t> rhs;
  real_t dt;
  virtual void forcings(real_t dt) = 0;

  arrvec_t<arr_2d_t> C;
  //interpolate from "skalar filed" velocity to  shifted grid courant field
  void intrp_courant()
  {
    arr_2d_t tmp, tmp_u;

    const rng_t &i = this->i;
    const rng_t &j = this->j;

    //filing halos for velocity filed
    this->xchng(u);
    this->xchng(w);
  
//    tmp_c = this->courant(0);
//    tmp = this->state(u);
    this->courant(0)(i+h,j) = .5*(this->state(u)(i^h,j)+this->state(u)((i^h)+1,j));
//    this->courant(0) = tmp_c;

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
   
      //interpolate velocity to fill courant field
      intrp_courant();

      //extrapolate courant field for mpdata (t+1/2)
      //(w tej kolejności mogę ciągle zaaplikowac prawą stronę do prędkości)

      forcings(dt / 2);
      this->xchng();
      this->advop(tht);
      this->advop(u); 
      this->advop(w);   
      this->cycle(); 
      forcings(dt / 2);
    }
  }
};

