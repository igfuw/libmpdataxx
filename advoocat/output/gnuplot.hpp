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

      static_assert(parent_t::n_dims < 3, "only 1D and 2D output supported");

      std::unique_ptr<Gnuplot> gp;

      void start(const int nt)
      {
        gp.reset(new Gnuplot());

        *gp 
	   << "set grid\n"
	   << "set border " << p.gnuplot_border << "\n"
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
	   << "set zrange " << p.gnuplot_zrange << "\n"
	   << "set xrange [0:" << this->mem->state(0).extent(0)-1 << "]\n"
	   << "set xlabel '" << p.gnuplot_xlabel << "'\n"
	   << "set ylabel '" << p.gnuplot_ylabel << "'\n"
	   << "set term svg dynamic\n"
        ;

        if (parent_t::n_dims == 1) // known at compile time
        {
          if (p.gnuplot_command == "splot") 
          {
            *gp << "set yrange [0:" << nt << "]\n";
            if (p.gnuplot_ylabel == "") *gp << "set ylabel 't/dt'\n";
          }
          
          *gp 
	     << "set output '" << p.gnuplot_output << "'\n"
             << p.gnuplot_command << " 0 notitle"
          ;

          for (int t = 0; t <= nt; t+=p.outfreq)
          {
	    for (const auto &v : p.outvars)
            {
	      *gp << ", '-'";
              if (p.gnuplot_command == "splot") *gp << " using 0:(" << t << "):1";
              *gp << " with " << p.gnuplot_with << " lt " << v.first << (
                t == 0 
                ? std::string(" title '") + v.second.name + "'"
                : std::string(" notitle")
              );
            }
          }
          *gp << "\n";
        }

        if (parent_t::n_dims == 2) // known at compile time
        {
          *gp 
	     << "set cbrange " << p.gnuplot_cbrange << "\n"
	     << (p.gnuplot_view != "map" ? "set pm3d at b\n" : "")
	     << "set yrange [0:" << this->mem->state(0).extent(1)-1 << "]\n"
	     << "set xtics out\n"
	     << "set ytics out\n"
	  ;
        }
      }

      void stop()
      {
        gp.reset();
      }
 
      // helper constructs to make it compilable for both 1D and 2D versions
      std::string binfmt(blitz::Array<typename parent_t::real_t, 1>) { assert(false); }
      std::string binfmt(blitz::Array<typename parent_t::real_t, 2> a) { return gp->binfmt(a); }

      void record(const int var)
      {
        if (parent_t::n_dims == 1) // known at compile time
        {
          gp->send(this->mem->state(var));
        }

        if (parent_t::n_dims == 2) // known at compile time
        {
	  *gp << "set output '" << boost::format(p.gnuplot_output) % this->outvars[var].name % this->n << "'\n";
	  *gp << "set title '"<< this->outvars[var].name << " @ t/dt=" << std::setprecision(3) << this->n << "'\n";
	  *gp << p.gnuplot_command << " '-' binary" << binfmt(this->mem->state(0)) 
	      << "with " << p.gnuplot_with << " notitle\n";
	  gp->sendBinary(this->mem->state(var).copy());
        }
      }

      public:

      struct params_t : parent_t::params_t 
      { 
	std::string 
          gnuplot_output,
          gnuplot_with = (
            parent_t::n_dims == 2 
	      ? std::string("image failsafe") // 2D
	      : std::string("lines")          // 1D
          ),
          gnuplot_command = std::string("splot"),
          gnuplot_xlabel = std::string("X"),
          gnuplot_ylabel = (
            parent_t::n_dims == 2 
              ? std::string("Y") // 2D
              : std::string("")  // 1D
          ),
          gnuplot_view = std::string(""), 
          gnuplot_zrange = std::string("[*:*]"),
          gnuplot_cbrange = std::string("[*:*]"),
          gnuplot_border = std::string(""); 
        int gnuplot_maxcolors = 100;
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
