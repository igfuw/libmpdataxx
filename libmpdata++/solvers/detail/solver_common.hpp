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

#include <libmpdata++/formulae/opts.hpp>

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

        enum { n_eqs = ct_params_t::n_eqs };
        enum { halo = minhalo }; 
        enum { n_dims = ct_params_t::n_dims };
        enum { n_tlev = n_tlev_ };

        typedef typename ct_params_t::real_t real_t;
        typedef blitz::Array<real_t, n_dims> arr_t;

	protected: 

        long long int timestep = 0;
        std::vector<int> n; 

        typedef concurr::detail::sharedmem<real_t, n_dims, n_tlev> mem_t; 
	mem_t *mem;

	// helper methods invoked by solve()
	virtual void advop(int e) = 0;

	virtual void cycle(int e) final
	{ 
	  n[e] = (n[e] + 1) % n_tlev - n_tlev;  // -n_tlev so that n+1 does not give out of bounds
          if (e == n_eqs - 1) this->mem->cycle(); 
	}

	virtual void xchng(int e, int l = 0) = 0; // TODO: make l -> -l

        private:
      
#if !defined(NDEBUG)
        bool 
          hook_ante_step_called = true, // initially true to handle nt=0 
          hook_post_step_called = true, // 
          hook_ante_loop_called = false, 
          hook_post_loop_called = false;
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

        virtual void hook_ante_loop(const int nt) 
        {
#if !defined(NDEBUG)
          hook_ante_loop_called = true;
          if (nt > 0)
          {
	    hook_ante_step_called = false;
	    hook_post_step_called = false;
          }
#endif
        }

// TODO: this conflicts with multiple solve() calls - consider removing it
        virtual void hook_post_loop() 
        {
#if !defined(NDEBUG)
          hook_post_loop_called = true;
#endif
        }

	public:

        struct rt_params_t 
        {
          std::array<int, n_dims> span;
        };

	// ctor
	solver_common(mem_t *mem, const rt_params_t &p) :
	  n(n_eqs, 0), 
          mem(mem)
	{
	  static_assert(n_eqs > 0, "!");
        }

        // dtor
        virtual ~solver_common()
        {
#if !defined(NDEBUG)
	  assert(hook_ante_step_called && "any overriding hook_ante_step() must call parent_t::hook_ante_step()");
	  assert(hook_post_step_called && "any overriding hook_post_step() must call parent_t::hook_post_step()");
	  assert(hook_ante_loop_called && "any overriding hook_ante_loop() must call parent_t::hook_ante_loop()");
	  assert(hook_post_loop_called && "any overriding hook_post_loop() must call parent_t::hook_post_loop()");
#endif
        }

	virtual void solve(int nt) final
	{   
          // multiple calls to sovlve() are meant to advance the solution by nt
          nt += timestep;

          // being generous about out-of-loop barriers 
          if (timestep == 0)
          {
	    this->mem->barrier();
	    hook_ante_loop(nt);
	    this->mem->barrier();
          }

	  while (timestep < nt)
	  {   
	    // progress-bar info through thread name (check top -H)
	    monitor(float(timestep) / nt); 

            // might be used to implement multi-threaded signal handling
            this->mem->barrier();
            if (this->mem->panic) break;

            // proper solver stuff
            hook_ante_step();
	    for (int e = 0; e < n_eqs; ++e) xchng(e);
	    for (int e = 0; e < n_eqs; ++e) advop(e);
	    for (int e = 0; e < n_eqs; ++e) cycle(e); // note: cycle assumes ascending loop index
            timestep++;
            hook_post_step();
	  }   

          this->mem->barrier();
          hook_post_loop();
          this->mem->barrier();
        }

        protected:

	// psi[n] getter - just to shorten the code
	virtual const arr_t &psi_n(int e) final
	{
	  return this->mem->psi[e][this->n[e]];
	}

        static rng_t rng_vctr(const int n) { return rng_t(0, n-1)^h^(halo-1); }
        static rng_t rng_sclr(const int n) { return rng_t(0, n-1)^halo; }
      };

      template<typename ct_params_t, int n_tlev, int minhalo, class enableif = void>
      class solver
      {}; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
