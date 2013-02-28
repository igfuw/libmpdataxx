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
    enum inhomo_e { euler, strang };

    template <class homo_solver, inhomo_e inhomo>
    class inhomo_solver : public homo_solver
    {
      public:

      typedef homo_solver parent_t;
      typedef typename parent_t::mem_t mem_t;

      protected:

      virtual void forcings(typename parent_t::real_t dt) = 0;

      typename parent_t::real_t dt;

      public:

      struct params_t : parent_t::params_t 
      { 
        typename parent_t::real_t dt; 
      };

      // 1D
      inhomo_solver(
	typename parent_t::mem_t *mem, 
	typename parent_t::bc_p &bcxl, 
	typename parent_t::bc_p &bcxr, 
	const rng_t &i, 
	const params_t &p
      ) :
	parent_t(mem, bcxl, bcxr, i, p), 
        dt(p.dt)
      {}

      // 2D
      inhomo_solver(
	typename parent_t::mem_t *mem, 
	typename parent_t::bc_p &bcxl, 
	typename parent_t::bc_p &bcxr, 
	typename parent_t::bc_p &bcyl, 
	typename parent_t::bc_p &bcyr, 
	const rng_t &i, 
	const rng_t &j, 
	const params_t &p
      ) :
	parent_t(mem, bcxl, bcxr, bcyl, bcyr, i, j, p), 
        dt(p.dt)
      {}

      // 3D
      inhomo_solver(
	typename parent_t::mem_t *mem, 
	typename parent_t::bc_p &bcxl, // TODO: encapsulating all these into a struct would shorten by 5 lines!
	typename parent_t::bc_p &bcxr, 
	typename parent_t::bc_p &bcyl, 
	typename parent_t::bc_p &bcyr, 
	typename parent_t::bc_p &bczl, 
	typename parent_t::bc_p &bczr, 
	const rng_t &i, // TODO: ditto
	const rng_t &j, 
	const rng_t &k, 
	const params_t &p
      ) :
	parent_t(mem, bcxl, bcxr, bcyl, bcyr, bczl, bczr, i, j, k, p), 
        dt(p.dt)
      {}

      void hook_ante_step()
      {
        switch (inhomo)
        {
          case euler: 
            forcings(dt); 
            break;
          case strang: 
            forcings(dt / 2); 
            break;
          default: 
            assert(false);
        }
      }

      void hook_post_step()
      {
        switch (inhomo)
        {
          case euler: 
            break;
          case strang: 
            forcings(dt / 2);
            break;
          default:
            assert(false);
        }
      } 
    };
  }; // namespace solvers
}; // namespace advoocat
