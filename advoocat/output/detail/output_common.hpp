/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

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

	private: 

	int n_out;

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
	    if (n % n_out == 0) record_all();
	  }
	  this->mem->barrier();
	}

	public:

	struct params_t : parent_t::params_t 
	{ 
	  int n_out; 
	  std::map<int, info> outvars;
	};

	// 2D ctor
	output_common(
	  typename parent_t::mem_t *mem,
	  typename parent_t::bc_p &bcxl,
	  typename parent_t::bc_p &bcxr,
	  typename parent_t::bc_p &bcyl,
	  typename parent_t::bc_p &bcyr,
	  const rng_t &i,
	  const rng_t &j,
	  const params_t &p
	) :
	parent_t(mem, bcxl, bcxr, bcyl, bcyr, i, j, p),
	  n_out(p.n_out), outvars(p.outvars)
	{}
      };
    }; // namespace detail
  }; // namespace output
}; // namespace advoocat
