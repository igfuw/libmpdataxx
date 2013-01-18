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
/*
- liczymy nowe psi
- ekstrapolujemy (łatwo, bo zawsze mamy zapisane stare psi)
- teraz pewnie czas wypełnić halo dla tego psi
- interpolowanie do brzegów oczek, żeby wyliczyć kuranty powinno 
  dać się zapisać jednym wyrażeniem blitzowym bez żadnych pętli, 
   bo prędkości potzebujemy tylko 
  "wewnątrz" obszaru "domena+halo" 
  (a po fill_halos w całym tym obszarze są sensowne prędkości)
*/

      //filing halos for velocity filed
      this->xchng(u);
      this->xchng(w);
    
      //interpolate to fill courant field

      //extrapolate courant field for mpdata
      // w tej kolejności mogę ciągle zaaplikowac prawą stronę do prędkości

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

