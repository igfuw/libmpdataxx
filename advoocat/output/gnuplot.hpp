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

#include <boost/format.hpp>

namespace advoocat
{
  namespace output
  {
    template <class solver_t>
    class gnuplot : public detail::output_timer<solver_t>
    {
      using parent_t = detail::output_timer<solver_t>;

      std::string binfmt;
      Gnuplot gp;

      void setup()
      {
        binfmt = gp.binfmt(this->mem->state(0));

        gp 
	   << "set grid\n"
	   << "set xlabel '" << p.gnuplot_xlabel << "'\n"
	   << "set ylabel '" << p.gnuplot_ylabel << "'\n"
           << "set xrange [0:" << this->mem->state(0).extent(0)-1 << "]\n"
           << "set yrange [0:" << this->mem->state(0).extent(1)-1 << "]\n"
          // progressive-rock connoisseur palette ;)
           << "set palette defined (0 '#ffffff', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n"
           << "set view " << p.gnuplot_view << "\n"
           << "set key font \",5\"\n "
           << "set contour base\n"
           << "set nosurface\n"
           << "set cntrparam levels 0\n";
      }
 
      void record(int var)
      {
        gp << "filename = '" << boost::format(p.gnuplot_output) % this->outvars[var].name % this->n << "'\n";
        gp << "set term svg dynamic name filename\n";
	gp << "set output filename\n";
        gp << "set title '"<< this->outvars[var].name << " @ t/dt=" << std::setprecision(3) << this->n << "'\n";
  //         << "set cbrange [298.5:302]\n"
        gp << "splot '-' binary" << binfmt << "with image notitle\n";
        gp.sendBinary(this->mem->state(var).copy());
      }

      public:

      struct params_t : parent_t::params_t 
      { 
	std::string 
          gnuplot_output,
          gnuplot_xlabel = "X",
          gnuplot_ylabel = "Y",
          gnuplot_view = "map"; 
      };

      const params_t p; // that's a copy - convenient but might be memory-consuming

      // ctor
      gnuplot(
	typename parent_t::ctor_args_t args,
	const params_t &p
      ) : parent_t(args, p), p(p)
      {}
    }; 
  }; // namespace output
}; // namespace advoocat
