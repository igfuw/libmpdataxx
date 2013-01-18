/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

template <typename real_t = float>
class cyclic_1d
{
  // member fields
  rng_t left_halo, rght_halo;
  rng_t left_edge, rght_edge;;

  public:

  // ctor
  cyclic_1d(
    const rng_t &i, int halo
  ) :
    left_halo(i.first()-halo,   i.first()-1    ),
    rght_edge(i.last() -halo+1, i.last()       ),
    rght_halo(i.last() +1,     i.last() +halo  ),
    left_edge(i.first(),       i.first()+halo-1)
  {} 

  // method invoked by the solver
  void fill_halos(const arr_1d_t<real_t> &a)
  {
    a(left_halo) = a(rght_edge);     
    a(rght_halo) = a(left_edge);     
  }
};

template<int d, typename real_t = float>
class cyclic_2d
{
  // member fields
  rng_t left_halo, rght_halo;
  rng_t left_edge, rght_edge;;

  public:

  // ctor
  cyclic_2d(
    const rng_t &i, const rng_t &j, int halo
  ) :
    left_halo(i.first()-halo, i.first()-1),
    rght_edge(i.last()-halo+1, i.last()  ),
    rght_halo(i.last()+1, i.last()+halo  ), 
    left_edge(i.first(), i.first()+halo-1)
  {} 

  // method invoked by the solver
  void fill_halos(const arr_2d_t<real_t> &a, const rng_t &j)
  {
    a(pi<d>(left_halo, j)) = a(pi<d>(rght_edge, j));     
    a(pi<d>(rght_halo, j)) = a(pi<d>(left_edge, j));     
  }
};

