/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/formulae/arakawa_c.hpp>
#include <libmpdata++/concurr/detail/sharedmem.hpp>

#include <libmpdata++/solvers/detail/monitor.hpp>

#include <libmpdata++/bcond/detail/bcond_common.hpp>

#include <array>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      using namespace libmpdataxx::arakawa_c;

      constexpr int max(const int a, const int b)
      {
        return a > b ? a : b;
      }

      template <typename ct_params_t, int n_tlev_, int minhalo>
      class solver_common
      {
	public:

        enum { n_eqns = ct_params_t::n_eqns };
        enum { halo = minhalo }; 
        enum { n_dims = ct_params_t::n_dims };
        enum { n_tlev = n_tlev_ };
        
        using ct_params_t_ = ct_params_t; // propagate ct_params_t mainly for output purposes
        typedef typename ct_params_t::real_t real_t;
        typedef blitz::Array<real_t, n_dims> arr_t;
        using bcp_t = std::unique_ptr<bcond::detail::bcond_common<real_t, halo>>;

        using ix = typename ct_params_t::ix;

        using advance_arg_t = typename std::conditional<ct_params_t::var_dt, real_t, int>::type;

	protected: 
        // TODO: output common doesnt know about ct_params_t
        bool var_dt = ct_params_t::var_dt;

        std::array<std::array<bcp_t, 2>, n_dims> bcs;

        const int rank;

        // di, dj, dk declared here for output purposes
        real_t prev_dt, dt, di, dj, dk, max_abs_div_eps, max_courant;
        std::array<real_t, n_dims> dijk;

	const idx_t<n_dims> ijk;

        long long int timestep = 0;
        real_t time = 0;
        std::vector<int> n; 

        typedef concurr::detail::sharedmem<real_t, n_dims, n_tlev> mem_t; 
	mem_t *mem;

	// helper methods invoked by solve()
	virtual void advop(int e) = 0;

	virtual void cycle(int e) final
	{ 
	  n[e] = (n[e] + 1) % n_tlev - n_tlev;  // -n_tlev so that n+1 does not give out of bounds
          if (e == n_eqns - 1) mem->cycle(rank); 
	}

	virtual void xchng(int e) = 0;
        // TODO: implement flagging of valid/invalid halo for optimisations

        virtual void xchng_vctr_alng(arrvec_t<arr_t>&, const bool ad = false) = 0;

        void set_bcs(const int &d, bcp_t &bcl, bcp_t &bcr)
        {
          bcs[d][0] = std::move(bcl);
          bcs[d][1] = std::move(bcr);
        }

	virtual real_t courant_number(const arrvec_t<arr_t>&) = 0;
	virtual real_t max_abs_vctr_div(const arrvec_t<arr_t>&) = 0;
       
        // return false if advector does not change in time
        virtual bool calc_gc() {return false;}

        virtual void scale_gc(const real_t time, const real_t cur_dt, const real_t old_dt) = 0;

        private:
      
#if !defined(NDEBUG)
        bool 
          hook_ante_step_called = true, // initially true to handle nt=0 
          hook_post_step_called = true, // 
          hook_ante_loop_called = true;
#endif

        protected:

        virtual void hook_ante_step() 
        { 
          // sanity check if all subclasses call their parents' hooks
#if !defined(NDEBUG)
          hook_ante_step_called = true;
#endif
        }

        virtual void hook_post_step() 
        {
#if !defined(NDEBUG)
          hook_post_step_called = true;
#endif
        }

        virtual void hook_ante_loop(const advance_arg_t nt) 
        {
#if !defined(NDEBUG)
          hook_ante_loop_called = true;
#endif
          // fill halos in velocity field
          this->xchng_vctr_alng(mem->GC);
         
          // adaptive timestepping - for constant in time velocity it suffices
          // to change the timestep once and do a simple scaling of advector
          if (ct_params_t::var_dt)
          {
            real_t cfl = courant_number(mem->GC);
            if (cfl > 0)
            {
              dt *= max_courant / cfl;
              scale_gc(time, dt, prev_dt);
            }
          }
        }

	public:

        const real_t time_() const { return time;}

        struct rt_params_t 
        {
          std::array<int, n_dims> grid_size;
          real_t dt=0, max_abs_div_eps = blitz::epsilon(real_t(44)), max_courant = 0.5;
        };

	// ctor
	solver_common(
          const int &rank, 
          mem_t *mem, 
          const rt_params_t &p, 
          const decltype(ijk) &ijk
        ) :
          rank(rank), 
          prev_dt(p.dt),
          dt(p.dt),
          di(0),
          dj(0),
          dk(0),
          max_abs_div_eps(p.max_abs_div_eps),
          max_courant(p.max_courant),
	  n(n_eqns, 0), 
          mem(mem),
          ijk(ijk)
	{
          // compile-time sanity checks
	  static_assert(n_eqns > 0, "!");

          // run-time sanity checks
          for (int d = 0; d < n_dims; ++d)
            if (p.grid_size[d] < 1) 
              throw std::runtime_error("bogus grid size");
        }

        // dtor
        virtual ~solver_common()
        {
#if !defined(NDEBUG)
	  assert(hook_ante_step_called && "any overriding hook_ante_step() must call parent_t::hook_ante_step()");
	  assert(hook_post_step_called && "any overriding hook_post_step() must call parent_t::hook_post_step()");
	  assert(hook_ante_loop_called && "any overriding hook_ante_loop() must call parent_t::hook_ante_loop()");
#endif
        }

	virtual void solve(advance_arg_t nt) final
	{   
          // multiple calls to sovlve() are meant to advance the solution by nt
          // TODO: does it really work with var_dt ? we do not advance by time exactly ...
          nt += ct_params_t::var_dt ? time : timestep;

          // being generous about out-of-loop barriers 
          if (timestep == 0)
          {
	    mem->barrier();
#if !defined(NDEBUG)
	    hook_ante_loop_called = false;
#endif
	    hook_ante_loop(nt);
	    mem->barrier();
          }

          // moved here so that if an exception is thrown from hook_ante_loop these do not cause complaints
#if !defined(NDEBUG)
	  hook_ante_step_called = false;
	  hook_post_step_called = false;
#endif
          // higher-order temporal interpolation for output requires doing a few additional steps
          int additional_steps = ct_params_t::out_intrp_ord;
	  while (ct_params_t::var_dt ? (time < nt || additional_steps > 0) : timestep < nt)
	  {   
	    // progress-bar info through thread name (check top -H)
	    monitor(float(ct_params_t::var_dt ? time : timestep) / nt);  // TODO: does this value make sanse with repeated advence() calls?

            // might be used to implement multi-threaded signal handling
            mem->barrier();
            if (mem->panic) break;

            // proper solver stuff
            
            // for variable in time velocity calculate advector at n+1/2, returns false if
            // velocity does not change in time
            bool var_gc = calc_gc();

            // for variable in time velocity with adaptive time-stepping modify advector
            // to keep the Courant number roughly constant
            if (var_gc && ct_params_t::var_dt)
            {
              real_t cfl = courant_number(mem->GC);
              if (cfl > 0)
              {
                do 
                {
                  dt *= max_courant / cfl;
                  calc_gc();
                  cfl = courant_number(mem->GC);
                }
                while (cfl > max_courant);
              }
            }
            
            hook_ante_step();

	    for (int e = 0; e < n_eqns; ++e) scale(e, ct_params_t::hint_scale(e));

	    for (int e = 0; e < n_eqns; ++e) xchng(e);
	    for (int e = 0; e < n_eqns; ++e) 
            { 
              advop(e);
              if (e != n_eqns - 1) mem->barrier();
            }
	    for (int e = 0; e < n_eqns; ++e) cycle(e); // note: cycle assumes ascending loop index

	    for (int e = 0; e < n_eqns; ++e) scale(e, -ct_params_t::hint_scale(e));

            timestep++;
            time = ct_params_t::var_dt ? time + dt : timestep * dt;
            prev_dt = dt;
            hook_post_step();

            if (time >= nt) additional_steps--;
	  }   

          mem->barrier();
          // note: hook_post_loop was removed as conficling with multiple-advance()-call logic
        }

        protected:

	// psi[n] getter - just to shorten the code
        // note that e.g. in hook_post_loop it points rather to 
        // psi^{n+1} than psi^{n} (hence not using the name psi_n)
	virtual arr_t &state(const int &e) final
	{
	  return mem->psi[e][n[e]];
	}

        static rng_t rng_vctr(const rng_t &rng) { return rng^h^(halo-1); }
        static rng_t rng_sclr(const rng_t &rng) { return rng^halo; }

        private:
        void scale(const int &e, const int &exp)
        {
          if (exp == 0) return;
          else if (exp > 0) state(e)(ijk) /= (1 << exp);
          else if (exp < 0) state(e)(ijk) *= (1 << -exp);
        }
      };

      template<typename ct_params_t, int n_tlev, int minhalo, class enableif = void>
      class solver
      {}; 
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
