/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "../blitz.hpp"
#include "../idxperm.hpp"
#include "../arakawa_c.hpp"


namespace advoocat
{
  namespace formulae
  {
    namespace mpdata 
    {
      using namespace arakawa_c;
      using idxperm::pi;
// TODO
/*
  struct varsgn
  {
    template<class T> auto aon(const T &x) -> decltype(abs(x))
    {
      return abs(x);
    }
  };

  struct posdef
  {
    template<class T> T aon(const T &x)
    {
      return x;
    }
  };
*/
    template<class nom_t, class den_t>
    inline auto frac(
      const nom_t &nom, const den_t &den
    ) return_macro(,
      where(den > 0, nom / den, 0)
    ) 

    // 1D
    template <class arr_1d_t>
    inline auto A(
      const arr_1d_t &psi, 
      const rng_t &i 
    ) return_macro(,
      frac(
	  abs(psi(i+1)) 
	- abs(psi(i  )),
	// ----------------------
	  abs(psi(i+1)) 
	+ abs(psi(i  ))
      ) 
    ) 

    // 2D
    template<int d, class arr_2d_t>
    inline auto A(
      const arr_2d_t &psi, 
      const rng_t &i, 
      const rng_t &j
    ) return_macro(,
      frac(
	  abs(psi(pi<d>(i+1, j))) 
	- abs(psi(pi<d>(i,   j))),
	// ----------------------
	  abs(psi(pi<d>(i+1, j))) 
	+ abs(psi(pi<d>(i,   j)))
      ) 
    ) 

    template<int d, class arr_2d_t>
    inline auto B(
      const arr_2d_t &psi, 
      const rng_t &i, 
      const rng_t &j
    ) return_macro(,
     frac(//aon=abs or not -> depending on the template paramater
	  abs(psi(pi<d>(i+1, j+1))) 
	+ abs(psi(pi<d>(i,   j+1))) 
	- abs(psi(pi<d>(i+1, j-1))) 
	- abs(psi(pi<d>(i,   j-1))),
	// ------------------------
	  abs(psi(pi<d>(i+1, j+1))) 
	+ abs(psi(pi<d>(i,   j+1))) 
	+ abs(psi(pi<d>(i+1, j-1))) 
	+ abs(psi(pi<d>(i,   j-1)))
      ) / 2
    )

    template<int d, class arr_2d_t>
    inline auto C_bar(
      const arr_2d_t &C, 
      const rng_t &i, 
      const rng_t &j
    ) return_macro(,
      (
	C(pi<d>(i+1, j+h)) + 
	C(pi<d>(i,   j+h)) +
	C(pi<d>(i+1, j-h)) + 
	C(pi<d>(i,   j-h)) 
      ) / 4
    )

    template<class arr_1d_t>
    inline auto antidiff(
      const arr_1d_t &psi, 
      const rng_t &i, 
      const arr_1d_t &C
    ) return_macro(,
      abs(C(i+h)) 
      * (1 - abs(C(i+h))) 
      * A(psi, i) 
    ) 

    template <int dim, class arr_2d_t>
    inline auto antidiff(
      const arr_2d_t &psi, 
      const rng_t &i, 
      const rng_t &j,
      const arrvec_t<arr_2d_t> &C
    ) return_macro(,
      abs(C[dim](pi<dim>(i+h, j))) 
      * (1 - abs(C[dim](pi<dim>(i+h, j)))) 
      * A<dim>(psi, i, j) 
      - C[dim](pi<dim>(i+h, j)) 
      * C_bar<dim>(C[dim-1], i, j)
      * B<dim>(psi, i, j)
    ) 

    };  //namespace mpdata
  };  //namespace formulae
}; // namespace advoocat

