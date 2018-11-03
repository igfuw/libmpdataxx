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
      trapez,  // paraphrase of trapezoidal rule:         psi^n+1 = ADV(psi^n + 1/2 * R^n) + 1/2 * R^n+1 
      mixed    // allows implementation of arbitrary mixed explicit/implicit schemes
    };

    const std::map<rhs_scheme_t, std::string> scheme2string {
      {euler_a, "euler_a"},
      {euler_b, "euler_b"},
      {trapez , "trapez" },
      {mixed  , "mixed"  }
    };
    
    struct mpdata_rhs_family_tag {};

    template <class ct_params_t, int minhalo = 0>
    class mpdata_rhs : public mpdata<ct_params_t, minhalo>
    {
      using parent_t = mpdata<ct_params_t, minhalo>;

      enum { n = 0 }; // just to make n, n+1 look nice :)

#if !defined(NDEBUG)
      bool update_rhs_called = true; // so that it nt=0 there's no complain
#endif

      protected:
      using solver_family = mpdata_rhs_family_tag;

      // member fields
      arrvec_t<typename parent_t::arr_t> &rhs;

      virtual void update_rhs(
        arrvec_t<typename parent_t::arr_t> &rhs, 
        const typename parent_t::real_t &dt,
        const int &at
      ) 
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

        // fill halos with data (e.g. for computing gradients)
        for (int e = 0; e < parent_t::n_eqns; ++e) this->xchng(e);
      }

      virtual void apply_rhs(
        const typename parent_t::real_t &dt_arg
      ) final
      {
        for (int e = 0; e < parent_t::n_eqns; ++e) 
        {
          // do nothing for equations with no rhs
          if (opts::isset(ct_params_t::hint_norhs, opts::bit(e))) continue;

          // otherwise apply the rhs
          this->state(e)(this->ijk) += dt_arg * rhs.at(e)(this->ijk);
        }
      }

      public:
      
      // ctor
      mpdata_rhs(
	typename parent_t::ctor_args_t args, 
	const typename parent_t::rt_params_t &p
      ) :
	parent_t(args, p), 
        rhs(args.mem->tmp[__FILE__][0])
      {
        assert(this->dt != 0);
      }

      // dtor
      ~mpdata_rhs()
      {
#if !defined(NDEBUG)
        assert(update_rhs_called && "any overriding update_rhs() must call parent_t::update_rhs()");
#endif
      }

      virtual void hook_mixed_rhs_ante_loop()
      {
        assert(false && "empty hook_mixed_rhs_ante_loop called");
      }

      virtual void hook_mixed_rhs_ante_step()
      {
        assert(false && "empty hook_mixed_rhs_ante_step called");
      }

      virtual void hook_mixed_rhs_post_step()
      {
        assert(false && "empty hook_mixed_rhs_post_step called");
      }

      void hook_ante_loop(const typename parent_t::advance_arg_t nt)
      {
        parent_t::hook_ante_loop(nt);

        switch ((rhs_scheme_t)ct_params_t::rhs_scheme)
        {
          case rhs_scheme_t::euler_a: 
          case rhs_scheme_t::euler_b: 
            break;
          case rhs_scheme_t::trapez:
            update_rhs(rhs, this->dt / 2, n);
            break;
          case rhs_scheme_t::mixed:
            hook_mixed_rhs_ante_loop();
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
            update_rhs(rhs, this->dt, n);
            break;
          case rhs_scheme_t::euler_b: 
            update_rhs(rhs, this->dt, n);
            apply_rhs(this->dt); 
            break;
          case rhs_scheme_t::trapez: 
            apply_rhs(this->dt / 2); 
            break;
          case rhs_scheme_t::mixed: 
            hook_mixed_rhs_ante_step();
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
            apply_rhs(this->dt);
            break;
          case rhs_scheme_t::euler_b: 
            break;
          case rhs_scheme_t::trapez: 
            update_rhs(rhs, this->dt / 2, n+1);
            apply_rhs(this->dt / 2);
            break;
          case rhs_scheme_t::mixed: 
            hook_mixed_rhs_post_step();
            break;
          default:
            assert(false);
        }
      } 

      static void alloc(
        typename parent_t::mem_t *mem, 
        const int &n_iters
      ) {
        // TODO: optimise to skip allocs for equations with no rhs
	parent_t::alloc(mem, n_iters);
        parent_t::alloc_tmp_sclr(mem, __FILE__, parent_t::n_eqns); // rhs array for each equation
      }
    };
  } // namespace solvers
} // namespace libmpdataxx
