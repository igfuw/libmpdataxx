/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "output_common.hpp"
#include <boost/timer/timer.hpp>

// TEMP!
#include <sstream>

namespace advoocat
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

	void hook_ante_loop()
	{
	  if (this->mem->rank() == 0) tmr.reset(new boost::timer::cpu_timer());
	  parent_t::hook_ante_loop();
	}

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
}; // namespace advoocat
