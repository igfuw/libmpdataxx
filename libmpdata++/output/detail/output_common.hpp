/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <map>

namespace libmpdataxx
{
  namespace output
  {
    namespace detail 
    {
      template <class solver_t>
      class output_common : public solver_t
      {
	using parent_t = solver_t;

        // TODO: not here!
        // sanity checks for ct_params_t
        static_assert(solver_t::n_dims > 0, "ct_params_t::n_dims missing?");
        static_assert(sizeof(typename solver_t::real_t), "ct_params_t::real_t missing?");
        static_assert(solver_t::n_eqs > 0, "ct_params_t::n_eqs missing?");

	protected:

	struct info { std::string name, unit; };
	std::map<int, info> outvars;

	int outfreq;

	virtual void record(const int var) {}
	virtual void start(const int nt) {}
	virtual void stop() {}

	void hook_ante_loop(const int nt)
	{
	  parent_t::hook_ante_loop(nt);
	  if (this->mem->rank() == 0) 
          {
            start(nt);
            record_all();
          }
	  this->mem->barrier();
	}

	void hook_post_loop()
	{
	  if (this->mem->rank() == 0) stop();
	  this->mem->barrier();
	  parent_t::hook_post_loop();
	}

	virtual void record_all()
	{
	  for (const auto &v : outvars) record(v.first);
	}

	void hook_post_step()
	{
	  parent_t::hook_post_step();
          // TODO: here the order or hooks matter -> not good :(

	  this->mem->barrier(); // waiting for all threads befor doing global output
	  if (this->mem->rank() == 0)
	  {
            // TODO: output of solver statistics every timesteps could probably go here
	    if (this->timestep % outfreq == 0) record_all();
	  }
	  this->mem->barrier(); // waiting for the output to be finished
	}

	public:

	struct rt_params_t : parent_t::rt_params_t 
	{ 
	  int outfreq = 1; 
	  std::map<int, info> outvars;
	};

	// ctor
	output_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
          parent_t(args, p),
	  outfreq(p.outfreq), 
          outvars(p.outvars)
	{
          // default value for outvars
          if (this->outvars.size() == 0 && parent_t::n_eqs == 1)
            outvars = {{0, {.name = "", .unit = ""}}};
        }
      };
    }; // namespace detail
  }; // namespace output
}; // namespace libmpdataxx
