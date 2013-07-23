/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

// TODO: move to solver::detail?

#include <libmpdata++/output/detail/output_common.hpp>
#include <boost/timer/timer.hpp>

// TEMP! TODO
#include <sstream>

namespace libmpdataxx
{
  namespace output
  {
    namespace detail 
    {
      template <class solver_t>
      class output_timer : public output_common<solver_t>
      {
	using parent_t = output_common<solver_t>;

        protected:

        std::unique_ptr<boost::timer::cpu_timer> tmr; 

        // TODO: move to ctor?
	void hook_ante_loop(const int nt)
	{
	  if (this->mem->rank() == 0) tmr.reset(new boost::timer::cpu_timer());
	  parent_t::hook_ante_loop(nt);
	}

        // TODO: move to dtor?
	void hook_post_loop()
	{
	  parent_t::hook_post_loop();
	  if (this->mem->rank() == 0) 
          {
            tmr->stop();
            boost::timer::cpu_times t = tmr->elapsed();
            std::ostringstream tmp;
            tmp << " wall time: " << double(t.wall) * 1e-9 << "s";
            tmp << " user time: " << double(t.user) * 1e-9 << "s";
            tmp << " system time: " << double(t.system) * 1e-9 << "s";
            std::cerr << tmp.str() << std::endl;
          }
	}

	public:

	// ctor
	output_timer(
	  typename parent_t::ctor_args_t args,
	  const typename parent_t::params_t &p
	) : parent_t(args, p)
	{}
      };
    }; // namespace detail
  }; // namespace output
}; // namespace libmpdataxx
