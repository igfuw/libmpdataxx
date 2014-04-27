/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief improved Euler inhomogeneous solver  
 *        (cf. eq. 32 in Smolarkiewicz 1998) // TODO cite
 */

#include <libmpdata++/solvers/mpdata.hpp>

#pragma once

namespace libmpdataxx
{
  namespace solvers
  {
    enum rhs_scheme_t 
    { 
      euler_a, // Euler's method, Eulerian spirit:        psi^n+1 = ADV(psi^n) + R^n
      euler_b, // Euler's method, semi-Lagrangian spirit: psi^n+1 = ADV(psi^n + R^n)
      trapez   // paraphrase of the Strang splitting:     psi^n+1 = ADV(psi^n + 1/2 * R^n) + 1/2 * R^n+1 
    };

    template <class ct_params_t>
    class mpdata_rhs : public mpdata<ct_params_t>
    {
      using parent_t = mpdata<ct_params_t>;

      enum { n = 0 }; // just to make n, n+1 look nice :)

#if !defined(NDEBUG)
      bool update_rhs_called = true; // so that it nt=0 there's no complain
#endif

      protected:

      // member fields
      typename parent_t::real_t dt;
      arrvec_t<typename parent_t::arr_t> &rhs;

//<listing-1>
      virtual void update_rhs(
        arrvec_t<typename parent_t::arr_t> &rhs, 
        const typename parent_t::real_t &dt,
        const int &at
      ) 
//</listing-1>
      {
        assert(at == n || at == n+1);
#if !defined(NDEBUG)
        update_rhs_called = true;
#endif
        // zero-out all rhs arrays
	for (int e = 0; e < parent_t::n_eqns; ++e) 
        {
          // do nothing for equations with no rhs
          if (opts::isset(ct_params_t::hint_norhs, opts::bit(e))) continue;

          // otherwise zero out the rhs
          rhs.at(e)(this->ijk) = 0;
        }
      }

      virtual void apply_rhs(
        const typename parent_t::real_t &dt
      ) final
      {
        for (int e = 0; e < parent_t::n_eqns; ++e) 
        {
          // do nothing for equations with no rhs
          if (opts::isset(ct_params_t::hint_norhs, opts::bit(e))) continue;

          // otherwise apply the rhs
          this->psi_n(e)(this->ijk) += dt * rhs.at(e)(this->ijk);
        }
      }

      public:

      struct rt_params_t : parent_t::rt_params_t 
      { 
        typename parent_t::real_t dt = 0; 
      };
      
      // ctor
      mpdata_rhs(
	typename parent_t::ctor_args_t args, 
	const rt_params_t &p
      ) :
	parent_t(args, p), 
        dt(p.dt),
        rhs(args.mem->tmp[__FILE__][0])
      {
        assert(dt != 0);
      }

      // dtor
      ~mpdata_rhs()
      {
#if !defined(NDEBUG)
        assert(update_rhs_called && "any overriding update_rhs() must call parent_t::update_rhs()");
#endif
      }

      void hook_ante_loop(int nt)
      {
        parent_t::hook_ante_loop(nt);

        switch ((rhs_scheme_t)ct_params_t::rhs_scheme)
        {
          case rhs_scheme_t::euler_a: 
          case rhs_scheme_t::euler_b: 
            break;
          case rhs_scheme_t::trapez:
            update_rhs(rhs, dt / 2, n);
            break;
          default: 
            assert(false);
        }
      }

      void hook_ante_step()
      {
        parent_t::hook_ante_step();

#if !defined(NDEBUG)
        update_rhs_called = false;
#endif

        switch ((rhs_scheme_t)ct_params_t::rhs_scheme)
        {
          case rhs_scheme_t::euler_a: 
            update_rhs(rhs, dt, n);
            break;
          case rhs_scheme_t::euler_b: 
            update_rhs(rhs, dt, n);
            apply_rhs(dt); 
            break;
          case rhs_scheme_t::trapez: 
            apply_rhs(dt / 2); 
            break;
          default: 
            assert(false);
        }
      }

      void hook_post_step()
      {
        parent_t::hook_post_step();
        switch ((rhs_scheme_t)ct_params_t::rhs_scheme)
        {
          case rhs_scheme_t::euler_a: 
            apply_rhs(dt);
            break;
          case rhs_scheme_t::euler_b: 
            break;
          case rhs_scheme_t::trapez: 
            update_rhs(rhs, dt / 2, n+1);
            apply_rhs(dt / 2);
            break;
          default:
            assert(false);
        }
      } 

      static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
      {
        // TODO: optimise to skip allocs for equations with no rhs
	parent_t::alloc(mem, p);
        parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, parent_t::n_eqns); // rhs array for each equation
      }
    };
  }; // namespace solvers
}; // namespace libmpdataxx
