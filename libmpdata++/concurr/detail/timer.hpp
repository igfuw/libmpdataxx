/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <boost/timer/timer.hpp>
#include <sstream>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail 
    {
      class timer 
      {
        std::unique_ptr<boost::timer::cpu_timer> tmr; 
        bool started = false;

        public:
      
        // ctor
        timer()
        {
          tmr.reset(new boost::timer::cpu_timer());
        }

        void resume()
        { 
          if (started) tmr->resume();
          else 
          {
            started = true;
            tmr->start();
          }
        }

        void stop()
        {
          tmr->stop();
        }

        void print()
        {
          boost::timer::cpu_times t = tmr->elapsed();
          std::ostringstream tmp;
          tmp << " wall time: " << double(t.wall) * 1e-9 << "s";
          tmp << " user time: " << double(t.user) * 1e-9 << "s";
          tmp << " system time: " << double(t.system) * 1e-9 << "s";
          std::cerr << tmp.str() << std::endl;
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
