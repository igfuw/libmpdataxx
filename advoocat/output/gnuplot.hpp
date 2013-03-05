/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "detail/output_common.hpp"

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

namespace advoocat
{
  namespace output
  {
    template <class solver_t>
    class gnuplot : public detail::output_common<solver_t>
    {
      using parent_t = detail::output_common<solver_t>;

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

      // 2D ctor
      gnuplot(
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
	plotfile(p.plotfile)
      {}
    }; 
  }; // namespace output
}; // namespace advoocat
