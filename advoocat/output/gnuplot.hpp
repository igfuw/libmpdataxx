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
      std::unique_ptr<Gnuplot> gp;

      void start()
      {
        gp.reset(new Gnuplot());
        binfmt = gp->binfmt(this->mem->state(0));

        *gp 
	   << "set grid\n"
	   << "set xlabel '" << p.gnuplot_xlabel << "'\n"
	   << "set ylabel '" << p.gnuplot_ylabel << "'\n"
           << "set xrange [0:" << this->mem->state(0).extent(0)-1 << "]\n"
           << "set yrange [0:" << this->mem->state(0).extent(1)-1 << "]\n"
           << "set zrange " << p.gnuplot_zrange << "\n"
           << "set cbrange " << p.gnuplot_cbrange << "\n"
           << "set xtics out\n"
           << "set ytics out\n"
           << "set palette defined ("
             "0 '#ffffff'," //         /\-
             "1 '#993399'," //        /  \-
             "2 '#00CCFF'," //  -----/    \---
             "3 '#66CC00'," // -----/      \---___
             "4 '#FFFF00'," //     /        \-    ---
             "5 '#FC8727'," //    /__________\-
             "6 '#FD0000'"  // 
           ") maxcolors " << p.gnuplot_maxcolors << "\n" 
           << "set view " << p.gnuplot_view << "\n"
           << "set border " << p.gnuplot_border << "\n"
           << "set key font \",5\"\n "
           << (p.gnuplot_view != "map" ? "set pm3d at b\n" : "")
           << "set cntrparam levels 0\n";
      }

      void stop()
      {
        gp.reset();
      }
 
      void record(int var)
      {
	*gp << "set output '" << boost::format(p.gnuplot_output) % this->outvars[var].name % this->n << "'\n";
        *gp << "set term svg dynamic\n";
        *gp << "set title '"<< this->outvars[var].name << " @ t/dt=" << std::setprecision(3) << this->n << "'\n";
  //         << "set cbrange [298.5:302]\n"
        *gp << "splot '-' binary" << binfmt << "with " << p.gnuplot_with << " notitle\n";
        gp->sendBinary(this->mem->state(var).copy());
      }

      public:

      struct params_t : parent_t::params_t 
      { 
	std::string 
          gnuplot_output,
          gnuplot_with = std::string("image failsafe"),
          gnuplot_xlabel = std::string("X"),
          gnuplot_ylabel = std::string("Y"),
          gnuplot_view = std::string(""), 
          gnuplot_zrange = std::string("[*:*]"),
          gnuplot_cbrange = std::string("[*:*]"),
          gnuplot_border = std::string(""); 
        int gnuplot_maxcolors = 100;
        //bool gnuplot_surface = false;
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
