/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "gnuplot-iostream_1d/test_gnuplot-iostream_1d.cpp"
 * \image html "../../tests/gnuplot-iostream_1d/figure.svg"
 */

#include <advoocat/lib.hpp>
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>
#include <tuple>

template <class T>
void setup(T &solver, int n, real_t C) {
  blitz::firstIndex i;
  solver.state() = exp(
    -sqr(i-n/2.) / (2.*pow(n/10, 2))
  );  
  solver.courant() = C; 
}

int main() {
  int n = 10, nt = 10;
  real_t C = .5;
  Gnuplot gp;

  gp << "set term svg size 1500,500 dynamic\n" 
     << "set output 'figure.svg'\n"     
     << "set multiplot layout 1,1\n" 
     << "set noxtics\n"
     << "set ztics .5\n"
     << "set xlabel 'X'\n"
     << "set ylabel 't'\n"
     << "set zrange [-1:1]\n"   
     << "set cbrange [-1:.5]\n"     
     << "set xyplane at -.5\n"
     << "set cbtics .5\n"
     << "set noytics\n"
     << "set palette maxcolors 18 defined (-1. \"red\", -.5 \"red\", -.5 \"blue\", -.17 \"green\", .16 \"brown\", .5 \"black\")\n";

  auto slvs = std::make_tuple(
    solvers::donorcell_1d<cyclic_1d>(n),
    solvers::mpdata_1d<2, cyclic_1d>(n),
    solvers::mpdata_1d<4, cyclic_1d>(n)
  );

  {
    auto solver = std::get<2>(slvs);
    setup(solver, n, C);

    // instructing gnuplot what to plot
    gp << "splot '-' using 0:(0):1 with lines notitle";
    for (int t = 0; t < nt; ++t) 
      gp << ", '-' using 0:(" << t+1 << "):1 with lines notitle";
    gp << "\n";

    // sending data for t=0
    gp.send(solver.state().copy());

    // timestepping and sending data for t>0
    for (int t = 0; t < nt; ++t)
    {
      solver.solve(1);
      gp.send(solver.state().copy());
    }
  } 




}
