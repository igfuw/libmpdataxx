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

	const typename parent_t::real_t outfreq;
        const std::string outdir;

	virtual void record(const int var) {}
	virtual void start(const typename parent_t::real_t) {}

	void hook_ante_loop(const typename parent_t::real_t tshift)
	{
	  parent_t::hook_ante_loop(tshift);
	  if (this->rank == 0) 
          {
            start(tshift);
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
            auto old_time = this->time - this->dt;
            // record if time and old_time fall into different outfreq length intervals 
	    if ((std::floor(this->time / outfreq) - std::floor(old_time / outfreq)) == 1) record_all();
          }
	  
	  this->mem->barrier(); // waiting for the output to be finished
	}

	public:

	struct rt_params_t : parent_t::rt_params_t 
	{ 
	  typename parent_t::real_t outfreq = 0; 
	  std::map<int, info_t> outvars;
          std::string outdir;
          // TODO: pass adiitional info? (command_line, library versions, ...)
	};

	// ctor
	output_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) :
          parent_t(args, p),
	  outfreq(p.outfreq), 
          outvars(p.outvars),
          outdir(p.outdir)
	{
          assert(outfreq >= this->dt);
          // default value for outvars
          if (this->outvars.size() == 0 && parent_t::n_eqns == 1)
            outvars = {{0, {.name = "", .unit = ""}}};
          
          // assign 1 to di, dj, dk for output purposes if they are not defined by the user
          for (auto ref : 
                std::vector<std::reference_wrapper<typename parent_t::real_t>>{this->di, this->dj, this->dk})
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
