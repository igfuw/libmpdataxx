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

        bool do_record = false;
	typename parent_t::real_t prev_time, record_time;
	const typename parent_t::advance_arg_t outfreq;
	const int outwindow;
        const std::string outdir;

        arrvec_t<typename parent_t::arr_t> &intrp_vars;

	virtual void record(const int var) {}
	virtual void start(const typename parent_t::advance_arg_t nt) {}
        
        typename parent_t::arr_t out_data(const int var)
        {
          return this->var_dt ? intrp_vars[var] : this->mem->advectee(var);
        }

	void hook_ante_loop(const typename parent_t::advance_arg_t nt)
	{
	  parent_t::hook_ante_loop(nt);

          if (this->var_dt)
          {
            for (const auto &v : outvars)
              intrp_vars[v.first](this->ijk) = this->mem->advectee(v.first)(this->ijk);
          }

	  if (this->rank == 0) 
          {
            record_time = this->time;
            start(nt);
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
         
          if (this->var_dt)
          {
            prev_time = this->time;
            auto next_time = this->time + this->dt;

            int next_idx = std::floor(next_time / outfreq);
            int curr_idx = std::floor(this->time / outfreq);

            if (next_idx > curr_idx) 
            {
              do_record = true;
              record_time = next_idx * outfreq;
              for (const auto &v : outvars)
                intrp_vars[v.first](this->ijk) = this->mem->advectee(v.first)(this->ijk);
            }
          }
        }

	void hook_post_step() 
	{
	  parent_t::hook_post_step();

	  this->mem->barrier(); // waiting for all threads befor doing global output

          if (this->var_dt && do_record)
          {
              for (const auto &v : outvars)
              {
                intrp_vars[v.first](this->ijk) *= (this->time - record_time) / (this->time - prev_time);
                intrp_vars[v.first](this->ijk) += (record_time - prev_time) / (this->time - prev_time) *
                                                   this->mem->advectee(v.first)(this->ijk);
              }
          }

	  if (this->rank == 0)
	  {
            //TODO: output of solver statistics every timesteps could probably go here

            if (do_record)
            {
              record_all();
              do_record = false;
            }
            //for (int t = 0; t < outwindow; ++t)
            //{
            //  if ((this->timestep - t) % outfreq == 0) record_all();
            //}
          }
	  
	  this->mem->barrier(); // waiting for the output to be finished
	}

	public:

	struct rt_params_t : parent_t::rt_params_t 
	{ 
	  typename parent_t::advance_arg_t outfreq = 1;
	  int outwindow = 1;
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
	  outwindow(p.outwindow),
          outvars(p.outvars),
          outdir(p.outdir),
          intrp_vars(args.mem->tmp[__FILE__][0])
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

        static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
        {
          parent_t::alloc(mem, n_iters);
          // TODO: only allocate for outvars !
          parent_t::alloc_tmp_sclr(mem, __FILE__, parent_t::n_eqns);
        }
      };
    } // namespace detail
  } // namespace output
} // namespace libmpdataxx
