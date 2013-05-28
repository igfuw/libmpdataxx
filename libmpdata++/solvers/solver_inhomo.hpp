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

namespace libmpdataxx
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

      private:

      bool update_forcings_called = false;

      protected:

      // member fields
      typename parent_t::real_t dt;

      virtual void update_forcings(
        arrvec_t<typename parent_t::arr_t> &rhs
      ) 
      {
#if !defined(NDEBUG)
        update_forcings_called = true;
#endif
        // zero-out all rhs arrays
	for (int e = 0; e < parent_t::n_eqs; ++e) 
          rhs.at(e)(this->ijk) = 0;
      }

      // this assumes explicit form - TODO: _explicit in a different file? (alloc, apply, update)
      virtual void apply_forcings(typename parent_t::real_t dt_arg)
      {
        for (int e = 0; e < parent_t::n_eqs; ++e) 
          this->state(e)(this->ijk) += dt_arg * rhs.at(e)(this->ijk);
      }

      arrvec_t<typename parent_t::arr_t> &rhs;

      public:

      struct params_t : parent_t::params_t 
      { 
        typename parent_t::real_t dt = 0; 
      };
      
      // ctor
      inhomo_solver(
	typename parent_t::ctor_args_t args, 
	const params_t &p
      ) :
	parent_t(args, p), 
        dt(p.dt),
        rhs(args.mem->tmp[__FILE__][0])
      {
        assert(dt != 0);
      }

      // dtor
      ~inhomo_solver()
      {
#if !defined(NDEBUG)
        assert(update_forcings_called && "any overriding update_forcings() must call parent_t::update_forcings()");
#endif
      }

      void hook_ante_loop(int nt)
      {
        parent_t::hook_ante_loop(nt);
        switch (inhomo)
        {
          case euler: 
            break;
          case strang:
            update_forcings(rhs);
            break;
          default: 
            assert(false);
        }
      }

      void hook_ante_step()
      {
        parent_t::hook_ante_step();
        switch (inhomo)
        {
          case euler: 
            update_forcings(rhs);
            apply_forcings(dt); 
            break;
          case strang: 
            apply_forcings(dt / 2); 
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
            update_forcings(rhs);
            apply_forcings(dt / 2);
            break;
          default:
            assert(false);
        }
      } 

      // TODO: merge the two allocs into one!

      // 1D version
      static void alloc(typename parent_t::mem_t *mem, const int nx)
      {
	parent_t::alloc(mem, nx);

        const rng_t i(0, nx-1);

	mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
	for (int e = 0; e < parent_t::n_eqs; ++e)
	  mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t(i^parent_t::halo)); 
      }

      // 2D version
      static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
      {
	parent_t::alloc(mem, nx, ny);

        const rng_t i(0, nx-1), j(0, ny-1);

	mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
	for (int e = 0; e < parent_t::n_eqs; ++e)
	  mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t(i^parent_t::halo, j^parent_t::halo));
      }
    };
  }; // namespace solvers
}; // namespace libmpdataxx
