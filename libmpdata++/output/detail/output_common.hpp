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

	typename parent_t::real_t prev_time;
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

	void hook_ante_step()
        {
          parent_t::hook_ante_step();
          prev_time = this->time;
        }

	void hook_post_step() 
	{
	  parent_t::hook_post_step();

	  this->mem->barrier(); // waiting for all threads befor doing global output
	  if (this->rank == 0)
	  {
            // TODO: output of solver statistics every timesteps could probably go here

            // assumes that we know what next dt will be !
            auto next_time = this->time + this->dt;

            // indices of outfreq length intervals into which next, current and previous timesteps fall
            int next_idx = std::floor(next_time / outfreq);
            int curr_idx = std::floor(this->time / outfreq);
            int prev_idx = std::floor(prev_time / outfreq);
            
            // to avoid repeated records
            bool record_all_called = false;
            
            // if next and current timestep fall into different intervals
            if (next_idx > curr_idx) 
            {
              auto next_diff = next_time - next_idx * outfreq;
              auto curr_diff = next_idx * outfreq - this->time;

              // check if current is closer to the interval boundary and if so record it
              if (curr_diff < next_diff)
              {
                record_all();
                record_all_called = true;
              }
            }

            // same logic as above - only for current and previous timestep, skip if record already done
            if (!record_all_called && curr_idx > prev_idx)
            {
              auto curr_diff = this->time - curr_idx * outfreq;
              auto prev_diff = curr_idx * outfreq - prev_time;
              
              if (curr_diff < prev_diff) record_all();
            }
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
