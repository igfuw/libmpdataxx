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

        typedef typename ct_params_t::real_t real_t;
        typedef blitz::Array<real_t, n_dims> arr_t;
        using bcp_t = std::unique_ptr<bcond::detail::bcond_common<real_t>>;

        using ix = typename ct_params_t::ix;

	protected: 
        
        const int rank;

        // declared here for output purposes
        real_t dt, di, dj, dk;

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

        virtual void hook_ante_loop(const real_t tshift) 
        {
#if !defined(NDEBUG)
          hook_ante_loop_called = true;
#endif
        }

	public:

        struct rt_params_t 
        {
          std::array<int, n_dims> grid_size;
          real_t dt=1;
        };

	// ctor
	solver_common(
          const int &rank, 
          mem_t *mem, 
          const rt_params_t &p, 
          const decltype(ijk) &ijk
        ) :
          rank(rank), 
          dt(p.dt),
          di(0),
          dj(0),
          dk(0),
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

	virtual void solve(real_t tshift) final
	{   
          // multiple calls to sovlve() are meant to advance the solution by nt
          tshift += time;

          // being generous about out-of-loop barriers 
          if (timestep == 0)
          {
	    mem->barrier();
#if !defined(NDEBUG)
	    hook_ante_loop_called = false;
#endif
	    hook_ante_loop(tshift);
	    mem->barrier();
          }

          // moved here so that if an exception is thrown from hook_ante_loop these do not cause complaints
#if !defined(NDEBUG)
	  hook_ante_step_called = false;
	  hook_post_step_called = false;
#endif

	  while (time < tshift)
	  {   
	    // progress-bar info through thread name (check top -H)
	    monitor(float(time) / tshift);  // TODO: does this value make sanse with repeated advence() calls?

            // might be used to implement multi-threaded signal handling
            mem->barrier();
            if (mem->panic) break;

            // proper solver stuff
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
            if (!ct_params_t::var_dt)
              time = timestep * dt;
            else
              time += dt;
            hook_post_step();
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

        static rng_t rng_vctr(const int &n) { return rng_t(0, n-1)^h^(halo-1); }
        static rng_t rng_sclr(const int &n) { return rng_t(0, n-1)^halo; }

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
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
