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

namespace advoocat
{
  namespace solvers
  {
    template <class homo_solver, bool naive = false>
    class inhomo_solver : public homo_solver
    {
      public:

      typedef homo_solver parent_t;
      typedef typename parent_t::mem_t mem_t;
      typedef typename mem_t::real_t real_t;

      protected:

      virtual void forcings(real_t dt) = 0;

      // psi getter
      typename parent_t::mem_t::arr_t psi(int e, int add = 0)
      {
	return this->mem.psi[e][this->mem.n[e] + add];
      }

      real_t dt;

      public:

      struct params_t : parent_t::params_t { real_t dt; };

      // 1D
      inhomo_solver(
	typename parent_t::mem_t &mem, 
	const rng_t &i, 
	const params_t &p
      ) :
	parent_t(mem, i, p), dt(p.dt)
      {}

      // 2D
      inhomo_solver(
	typename parent_t::mem_t &mem, 
	const rng_t &i, 
	const rng_t &j, 
	const params_t &p
      ) :
	parent_t(mem, i, j, p), dt(p.dt)
      {}

      // 3D
      inhomo_solver(
	typename parent_t::mem_t &mem, 
	const rng_t &i, 
	const rng_t &j, 
	const rng_t &k, 
	const params_t &p
      ) :
	parent_t(mem, i, j, k, p), dt(p.dt)
      {}

      void solve(int nt)
      {
	for (int t = 0; t < nt; ++t)
	{
	  if (!naive) forcings(dt / 2);
	  else forcings(dt);

	  parent_t::solve(1);

	  if (!naive) forcings(dt / 2);
	}
      }
    };

    template <class homo_solver>
    using inhomo_solver_naive = inhomo_solver<homo_solver, true>;
  }; // namespace solvers
}; // namespace advoocat
