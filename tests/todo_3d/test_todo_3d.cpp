/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// (<> should be used instead of "" in normal usage)
//#include "advoocat/mpdata_3d.hpp"
#include "advoocat/donorcell_3d.hpp"
#include "advoocat/cyclic_3d.hpp"
#include "advoocat/equip.hpp"

enum {x, y, z};

int main() 
{
  int n[] = {24, 24, 24};
  {
    equip<
      solvers::donorcell_3d<
        cyclic_3d<x>, 
        cyclic_3d<y>,
        cyclic_3d<z>,
        sharedmem_3d<>
      >
    > slv(n[x], n[y], n[z]);
  } 
/*
  {
    solvers::mpdata_3d<
      2,
      cyclic_3d<x>, 
      cyclic_3d<y>,
      cyclic_3d<z>
    > solver(n[x], n[y], n[z]);
  } 
  {
    solvers::mpdata_3d<
      44,
      cyclic_3d<x>, 
      cyclic_3d<y>,
      cyclic_3d<z>
    > solver(n[x], n[y], n[z]);
  } 
*/
}
