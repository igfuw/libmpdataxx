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

// TODO: hook_ante_loop() - zero rhs arrays

      void hook_ante_step()
      {
        parent_t::hook_ante_step();
        switch (inhomo)
        {
          case euler: 
//            apply_forcings(dt); 
            break;
          case strang: 
//            apply_forcings(dt / 2); 
            break;
          default: 
            assert(false);
        }
      }

      void hook_post_step()
      {
        parent_t::hook_post_step();
        switch (inhomo)
        {
          case euler: 
            break;
          case strang: 
//            update_forcings();
//            apply_forcings(dt / 2);
            break;
          default:
            assert(false);
        }
      } 

      
// TODO: member field with rhs arrays (arrvec_t)
// TODO: ctor (a in mpdata solvers)

      static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
      {
	parent_t::alloc(mem, nx, ny);

        // (i,j)-sized temporary fields
	mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());

	const rng_t i(0, nx-1), j(0, ny-1);
	for (int n=0; n < parent_t::n_eqs; ++n)
	  mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t( i, j ));
      }
    };
  }; // namespace solvers
}; // namespace advoocat
