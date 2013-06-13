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

	protected:

	int n =0;
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

	void record_all()
	{
	  for (const auto &v : outvars) record(v.first);
	}

	void hook_post_step()
	{
	  parent_t::hook_post_step();

	  if (this->mem->rank() == 0)
	  {
	    n++;
// TODO: output of solver statistics every timesteps could probably go here
	    if (n % outfreq == 0) record_all();
	  }
	  this->mem->barrier();
	}

	public:

	struct params_t : parent_t::params_t 
	{ 
	  int outfreq = 1; 
	  std::map<int, info> outvars;
	};

	// ctor
	output_common(
	  typename parent_t::ctor_args_t args,
	  const params_t &p
	) :
	parent_t(args, p),
	  outfreq(p.outfreq), outvars(p.outvars)
	{}
      };
    }; // namespace detail
  }; // namespace output
}; // namespace libmpdataxx
