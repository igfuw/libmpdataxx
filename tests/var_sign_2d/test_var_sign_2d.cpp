/** 
 * @file
 * @example var_sign_2d/test_var_sign_2d.cpp
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief test for advection of variable sign field 
 *
 * \include "var_sign_2d/test_var_sign_2d.cpp"
 * \image html "../../tests/var_sign_2d/figure.svg"
 */

// (<> should be used instead of "" in normal usage)
#include "advoocat/solvers/donorcell_2d.hpp"
#include "advoocat/solvers/mpdata_2d.hpp"
#include "advoocat/bcond/cyclic_2d.hpp"
#include "advoocat/equip.hpp"

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>
#include <vector>
enum {x, y};

template <class T>
void setup(T &solver, int n[2], typename T::real_t offset) 
{
  blitz::firstIndex i;
  blitz::secondIndex j;
  solver.state() = offset + exp(
    -sqr(i-n[x]/2.) / (2.*pow(n[x]/10, 2))
    -sqr(j-n[y]/2.) / (2.*pow(n[y]/10, 2))
  );  
  solver.courant(x) = .5; 
  solver.courant(y) = .25;
}

int main() 
{
  using namespace advoocat;

  int n[] = {24, 24}, nt = 96;
  Gnuplot gp;
  gp << "set term svg size 600,1100 dynamic fsize 5\n" 
     << "set output 'figure.svg'\n"     
     << "set multiplot layout 4,2 columnsfirst\n" 
     << "set border 4095\n"
     << "set xtics out\n"
     << "set ytics out\n"
     << "set xlabel 'X'\n"
     << "set ylabel 'Y'\n"
     << "set xrange [0:" << n[x]-1 << "]\n"   
     << "set yrange [0:" << n[y]-1 << "]\n"   
     << "set zrange [-.666:1]\n"   
     << "set palette maxcolors 42\n"
     << "set pm3d at b\n";

  std::string binfmt;
  using mem_t = sharedmem_2d<>;

  for (auto &offset : std::vector<float>({0,-.5})){
    gp << "set cbrange ["<< offset -.025<<":"<<offset+1.025<<"]\n";
    {
      equip<
        solvers::donorcell_2d<bcond::cyclic_2d<x>, bcond::cyclic_2d<y>, mem_t>
      > slv(n[x],n[y]);

      setup(slv, n, offset);
      binfmt = gp.binfmt(slv.state());
      gp << "set title 't=0'\n"
         << "splot '-' binary" << binfmt
         << "with lines notitle\n";
      gp.sendBinary(slv.state().copy());
      slv.advance(nt);
      gp << "set title 'donorcell t="<<nt<<"'\n"
         << "splot '-' binary" << binfmt
         << "with lines notitle\n";
      gp.sendBinary(slv.state().copy());
    } {
      const int it = 2;
      equip<
        solvers::mpdata_2d<it, bcond::cyclic_2d<x>, bcond::cyclic_2d<y>, mem_t>
      > slv(n[x],n[y]); 
      setup(slv, n, offset); 
      slv.advance(nt);
      gp << "set title 'mpdata<" << it << "> "
         << "t=" << nt << "'\n"
         << "splot '-' binary" << binfmt
         << "with lines notitle\n";
      gp.sendBinary(slv.state().copy());
    } {
      const int it = 4;
      equip<
        solvers::mpdata_2d<it, bcond::cyclic_2d<x>, bcond::cyclic_2d<y>, mem_t>
      > slv(n[x],n[y]); 
      setup(slv, n, offset); 
      slv.advance(nt); 
      gp << "set title 'mpdata<" << it << "> "
         << "t=" << nt << "'\n"
         << "splot '-' binary" << binfmt
         << "with lines notitle\n";
      gp.sendBinary(slv.state().copy());
    }
  }
}
