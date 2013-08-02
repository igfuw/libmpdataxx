/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/arakawa_c.hpp>
#include <libmpdata++/concurr/detail/sharedmem.hpp>

#include <libmpdata++/solvers/adv/detail/monitor.hpp>

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

      template <typename real_t_, int n_dims_, int n_eqs_, int n_tlev_, int minhalo>
      class solver_common
      {
	public:

        // using enums as "public static const int" would need instantiation
        enum { halo = minhalo }; 
        enum { n_dims = n_dims_ };
        enum { n_eqs = n_eqs_ }; // TODO: this does not need to be a template parameter!
        enum { n_tlev = n_tlev_ };

        typedef real_t_ real_t;
        typedef blitz::Array<real_t_, n_dims_> arr_t;

	void cycle_all()
	{ 
	  for (int e = 0; e < n_eqs; ++e) cycle(e);
	}


	protected: 

        long long int timestep = 0;
        std::vector<int> n; // TODO: why not std::array?

        typedef concurr::detail::sharedmem<real_t, n_dims, n_eqs, n_tlev> mem_t;
	mem_t *mem;

	// helper methods invoked by solve()
	virtual void advop(int e) = 0;
	void advop_all()
	{
	  for (int e = 0; e < n_eqs; ++e) advop(e);
	}

	void cycle(int e) 
	{ 
	  n[e] = (n[e] + 1) % n_tlev - n_tlev;  // -n_tlev so that n+1 does not give out of bounds
          if (e == n_eqs - 1) this->mem->cycle(); 
	}

	virtual void xchng(int e, int l = 0) = 0; // TODO: make l -> -l
	void xchng_all() 
	{   
	  for (int e = 0; e < n_eqs; ++e) xchng(e);
	}

        private:
      
#if !defined(NDEBUG)
        bool 
          hook_ante_step_called = true, // initially true to handle nt=0 
          hook_post_step_called = true, // 
          hook_ante_loop_called = false, 
          hook_post_loop_called = false;
#endif

	public:

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

        virtual void hook_post_loop() 
        {
#if !defined(NDEBUG)
          hook_post_loop_called = true;
#endif
        }

	// ctor
	solver_common(mem_t *mem) :
	  n(n_eqs, 0), 
          mem(mem)
	{}

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

	virtual void solve(const int nt) final
	{   
          // being generous about out-of-loop barriers 
          this->mem->barrier();
          hook_ante_loop(nt);
          this->mem->barrier();

	  while (timestep < nt)
	  {   
	    // progress-bar info through thread name (check top -H)
	    monitor(float(timestep) / nt); 

            // multi-threaded SIGTERM handling
            this->mem->barrier();
            if (this->mem->panic) break;

            // proper solver stuff
            hook_ante_step();
	    xchng_all();
	    advop_all();
	    cycle_all();
            timestep++;
            hook_post_step();
	  }   

          this->mem->barrier();
          hook_post_loop();
          this->mem->barrier();
        }

	// psi getter
	arr_t state(int e, int add = 0) // TODO: rename it to something like psi() or subdomain_with_halos()!
	{
	  return this->mem->psi[e][this->n[e] + add];
	}

        protected:

        static rng_t rng_vctr(const int n) { return rng_t(0, n-1)^h^(halo-1); }
        static rng_t rng_sclr(const int n) { return rng_t(0, n-1)^halo; }
      };

      template<typename real_t, int n_dims, int n_eqs, int n_tlev, int minhalo>
      class solver
      {}; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
