/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief a thin wrapper around the gnuplot-iostream library intended 
 *   for debugging/demonstration purposes only (used in many tests)
 */

#pragma once

#include <libmpdata++/output/detail/output_timer.hpp>

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

#include <boost/format.hpp>

namespace libmpdataxx
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

        // some common 1D/2D settings
        *gp 
	   << "set grid\n"
	   << "set border " << p.gnuplot_border << "\n"
	   << "set palette " /*defined (" // makes gnuplot discard maxcolors :(
	     "0 '#ffffff'," //         /\-
	     "1 '#993399'," //        /  \-
	     "2 '#00CCFF'," //  -----/    \---
	     "3 '#66CC00'," // -----/      \---___
	     "4 '#FFFF00'," //     /        \-    ---
	     "5 '#FC8727'," //    /__________\-
	     "6 '#FD0000'"  // 
	   ")*/ << " maxcolors " << p.gnuplot_maxcolors << "\n" 
	   << "set view " << p.gnuplot_view << "\n"
	   << "set zrange " << p.gnuplot_zrange << "\n"
	   << "set xrange [0:" << this->mem->state(0).extent(0) << "]\n"
	   << "set xlabel '" << p.gnuplot_xlabel << "'\n"
	   << "set ylabel '" << p.gnuplot_ylabel << "'\n"
	   << "set term " << p.gnuplot_term << "\n"
        ;

        // 1D settings
        if (parent_t::n_dims == 1) // known at compile time
        {
          if (p.gnuplot_command == "splot") 
          {
            *gp << "set yrange [0:" << nt << "]\n";
            assert(p.gnuplot_yrange == "[*:*]" && "gnupot_yrange was specified for a 1D splot where Y axis represents time");

            if (p.gnuplot_ylabel == "") *gp << "set ylabel 't/dt'\n";
          }
          else if (p.gnuplot_command == "plot") 
          {
            *gp << "set yrange " << p.gnuplot_yrange << "\n";
          } 
          else assert(false);
          
          *gp 
	     << "set output '" << p.gnuplot_output << "'\n"
             << p.gnuplot_command << " 1/0 notitle"
          ;

          for (int t = 0; t <= nt; t+=p.outfreq)
          {
	    for (const auto &v : p.outvars)
            {
	      *gp << ", '-'";
              if (p.gnuplot_command == "splot") *gp << " using 0:(" << t << "):1";
              *gp << " with " << p.gnuplot_with;
 
              *gp << " lt ";
              if (p.outvars.size() == 1) *gp <<  p.gnuplot_lt;
              else *gp << v.first;

              *gp << (
                t == 0 
                ? std::string(" title '") + v.second.name + "'"
                : std::string(" notitle")
              );
            }
          }
          *gp << "\n";
        }

        // 2D settings
        if (parent_t::n_dims == 2) // known at compile time
        {
          *gp 
	     << "set cbrange " << p.gnuplot_cbrange << "\n"
	     << "set yrange [0:" << this->mem->state(0).extent(1) << "]\n"
	     << "set xtics out\n"
	     << "set ytics out\n"
	     << "set size square\n"
	     << (p.gnuplot_surface ? "set" : "unset") << " surface\n"
	     << (p.gnuplot_contour ? "set" : "unset") << " contour\n"
	  ;
          assert(p.gnuplot_yrange == "[*:*]" && "gnupot_yrange was specified for a 2D splot where Y axis represents spatial dimension");

          if (p.gnuplot_contour)
          {
            *gp 
               << "unset clabel\n"
               << "set cntrparam " << p.gnuplot_cntrparam << "\n";
          }
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
	  *gp << p.gnuplot_command;
          bool imagebg = (p.gnuplot_with == "lines");
          if (imagebg)
          {
            float zmin, zmax;
            int count = sscanf(p.gnuplot_zrange.c_str(), "[%g:%g]", &zmin, &zmax);
            if (count != 2) zmin = 0;
            *gp << " '-' binary " << binfmt(this->mem->state(0))
                << " origin=(.5,.5," << zmin << ")" // TODO: dx/2, dy/2, 
	        << " with image failsafe notitle,";
          }
          *gp << " '-'" 
              << " binary" << binfmt(this->mem->state(0)) 
              << " origin=(.5,.5,0)" // TODO: dx/2, dy/2, ?
	      << " with " << p.gnuplot_with << " lt " << p.gnuplot_lt << " notitle\n";
	  gp->sendBinary(this->mem->state(var).copy());
          if (imagebg) gp->sendBinary(this->mem->state(var).copy());
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
	      : std::string("lines ")         // 1D // TODO: histogram steps would be better than lines
          ),
          gnuplot_command = std::string("splot"),
          gnuplot_xlabel = std::string("x/dx"),
          gnuplot_ylabel = (
            parent_t::n_dims == 2 
              ? std::string("y/dy") // 2D
              : std::string("")  // 1D
          ),
          gnuplot_view = std::string(""), 
          gnuplot_zrange = std::string("[*:*]"),
          gnuplot_yrange = std::string("[*:*]"),
          gnuplot_cbrange = std::string("[*:*]"),
          gnuplot_border = std::string(""),
          gnuplot_lt = std::string("-1"), // black
          gnuplot_cntrparam = std::string(""),
          gnuplot_term = std::string("svg dynamic");
        int gnuplot_maxcolors = 100; 
        bool 
          gnuplot_contour = false,
          gnuplot_surface = true;
      };

      const params_t p; // TODO: that's a copy - convenient but might be memory-consuming, make a struct p.gnupot that would be copied

      // ctor
      gnuplot(
	typename parent_t::ctor_args_t args,
	const params_t &p
      ) : parent_t(args, p), p(p)
      {}
    }; 
  }; // namespace output
}; // namespace libmpdataxx
