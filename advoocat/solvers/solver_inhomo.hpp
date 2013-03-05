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
	typename parent_t::ctor_args_t args, 
	const params_t &p
      ) :
	parent_t(args, p), 
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
