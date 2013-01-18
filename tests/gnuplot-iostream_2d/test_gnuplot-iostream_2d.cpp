/** 
 * @file
 * @example gnuplot-iostream_2d/test_gnuplot-iostream_2d.cpp
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "gnuplot-iostream_2d/test_gnuplot-iostream_2d.cpp"
 * \image html "../../tests/gnuplot-iostream_2d/figure.svg"
 */

// (<> should be used instead of "" in normal usage)
#include "advoocat/mpdata_2d.hpp"
#include "advoocat/cyclic.hpp"

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>
enum {x, y};

template <class T>
void setup(T &solver, int n[2]) {
  blitz::firstIndex i;
  blitz::secondIndex j;
  solver.state() = exp(
    -sqr(i-n[x]/2.) / (2.*pow(n[x]/10, 2))
    -sqr(j-n[y]/2.) / (2.*pow(n[y]/10, 2))
  );  
  solver.courant(x) = .5; 
  solver.courant(y) = .25;
}

int main() {
  int n[] = {24, 24}, nt = 96;
  Gnuplot gp;
  gp << "set term svg size 500,1500 dynamic\n" 
     << "set output 'figure.svg'\n"     
     << "set multiplot layout 4,1\n" 
     << "set border 4095\n"
     << "set xtics out\n"
     << "set ytics out\n"
     << "unset ztics\n"    
     << "set xlabel 'X'\n"
     << "set ylabel 'Y'\n"
     << "set xrange [0:" << n[x]-1 << "]\n"   
     << "set yrange [0:" << n[y]-1 << "]\n"   
     << "set zrange [-.666:1]\n"   
     << "set cbrange [-.025:1.025]\n"     
     << "set palette maxcolors 42\n"
     << "set pm3d at b\n";
  std::string binfmt;
  {
    solvers::donorcell_2d<cyclic_2d<x>, cyclic_2d<y>> solver(n[x], n[y]);
    setup(solver, n);
    binfmt = gp.binfmt(solver.state());
    gp << "set title 't=0'\n"
       << "splot '-' binary" << binfmt
       << "with lines notitle\n";
    gp.sendBinary(solver.state().copy());
    solver.solve(nt);
    gp << "set title 'donorcell t="<<nt<<"'\n"
       << "splot '-' binary" << binfmt
       << "with lines notitle\n";
    gp.sendBinary(solver.state().copy());
  } 
  {
    const int it = 2;
    solvers::mpdata_2d<it, cyclic_2d<x>, cyclic_2d<y>> solver(n[x], n[y]); 
    setup(solver, n); 
    solver.solve(nt);
    gp << "set title 'mpdata<" << it << "> "
       << "t=" << nt << "'\n"
       << "splot '-' binary" << binfmt
       << "with lines notitle\n";
    gp.sendBinary(solver.state().copy());
  } 
  {
    const int it = 4;
    solvers::mpdata_2d<it, cyclic_2d<x>, cyclic_2d<y>> solver(n[x], n[y]); 
    setup(solver, n); 
    solver.solve(nt); 
    gp << "set title 'mpdata<" << it << "> "
       << "t=" << nt << "'\n"
       << "splot '-' binary" << binfmt
       << "with lines notitle\n";
    gp.sendBinary(solver.state().copy());
  }
}
