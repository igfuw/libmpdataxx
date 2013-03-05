/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "detail/output_timer.hpp"

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

namespace advoocat
{
  namespace output
  {
    template <class solver_t>
    class gnuplot : public detail::output_timer<solver_t>
    {
      using parent_t = detail::output_timer<solver_t>;

      std::string plotfile;

      void record(int var)
      {
std::cerr << "aqq " << this->n << " var=" << var << std::endl;
      }

      public:

      struct params_t : parent_t::params_t 
      { 
	std::string plotfile; 
      };

      // ctor
      gnuplot(
	typename parent_t::ctor_args_t args,
	const params_t &p
      ) :
      parent_t(args, p),
	plotfile(p.plotfile)
      {}
    }; 
  }; // namespace output
}; // namespace advoocat
