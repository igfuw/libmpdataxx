/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

template <class homo_solver>
class pressure_solver : public homo_solver
{
  enum{tht, prs, u, w};

  using real_t = typename homo_solver::real_t;

  //arrvec_t<typename homo_solver::arr_t> rhs;
  real_t dt;
  virtual void forcings(real_t dt) = 0;

  public:
  // ctor
  pressure_solver(int nx, int ny, real_t dt) :
    homo_solver(nx, ny), dt(dt)
  {}

  void solve(int nt)
  {
    for (int t = 0; t < nt; ++t)
    {
      //filing halos for velocity filed
      this->xchng(u);
      this->xchng(w);
    
      //interpolate velocity to fill courant field


      //extrapolate courant field for mpdata (t+1/2)
      //(w tej kolejności mogę ciągle zaaplikowac prawą stronę do prędkości)


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

