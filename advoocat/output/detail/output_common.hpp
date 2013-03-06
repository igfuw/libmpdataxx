/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <map>

namespace advoocat
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

	virtual void record(int var) {}
	virtual void setup() {}

	void hook_ante_loop()
	{
	  if (this->mem->rank() == 0) setup();
	  this->mem->barrier();
	  parent_t::hook_ante_loop();
	}

	void record_all()
	{
	  for (const auto &v : outvars) record(v.first);
	}

	void hook_ante_step()
	{
	  parent_t::hook_ante_step();
	  if (this->mem->rank() == 0)
	  {
	    if (n == 0) record_all();
	  }
	  this->mem->barrier();
	}

	void hook_post_step()
	{
	  parent_t::hook_post_step();
	  if (this->mem->rank() == 0)
	  {
	    n++;
	    if (n % outfreq == 0) record_all();
	  }
	  this->mem->barrier();
	}

	public:

	struct params_t : parent_t::params_t 
	{ 
	  int outfreq; 
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
}; // namespace advoocat
