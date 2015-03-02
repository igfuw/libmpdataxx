/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <map>
#include <vector>
#include <functional>

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

	struct info_t { std::string name, unit; };
	std::map<int, info_t> outvars;

	int outfreq;
	int outwindow;

	virtual void record(const int var) {}
	virtual void start(const int nt) {}

	void hook_ante_loop(const int nt) 
	{
	  parent_t::hook_ante_loop(nt);
	  if (this->rank == 0) 
          {
            start(nt);
            record_all();
          }
	  this->mem->barrier();
	}

	virtual void record_all()
	{
	  for (const auto &v : outvars) record(v.first);
	}

	void hook_post_step() 
	{
	  parent_t::hook_post_step();

	  this->mem->barrier(); // waiting for all threads befor doing global output
	  if (this->rank == 0)
	  {
            // TODO: output of solver statistics every timesteps could probably go here
            for (int t = 0; t < outwindow; ++t)
            {
	      if ((this->timestep - t) % outfreq == 0) record_all();
            }
          }
	  
	  this->mem->barrier(); // waiting for the output to be finished
	}

	public:

	struct rt_params_t : parent_t::rt_params_t 
	{ 
	  int outfreq = 1; 
	  int outwindow = 1;
	  std::map<int, info_t> outvars;
	};

	// ctor
	output_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
          parent_t(args, p),
	  outfreq(p.outfreq), 
	  outwindow(p.outwindow),
          outvars(p.outvars)
	{
          // default value for outvars
          if (this->outvars.size() == 0 && parent_t::n_eqns == 1)
            outvars = {{0, {.name = "", .unit = ""}}};
          
          // assign 1 to dt, di, dj, dk for output purposes if they are not defined by the user
          for (auto ref : 
                std::vector<std::reference_wrapper<typename parent_t::real_t>>{this->dt, this->di, this->dj, this->dk})
          {
            if (ref.get() == 0)
            {
              ref.get() = 1;
            }
          }
        }
      };
    }; // namespace detail
  }; // namespace output
}; // namespace libmpdataxx
