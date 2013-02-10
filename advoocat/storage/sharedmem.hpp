/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "../blitz.hpp"
#include <vector>
#include <unordered_map>

namespace advoocat
{
  template <int n_dims, int n_eqs_, typename real_t_>
  struct sharedmem_common
  {
    static_assert(n_eqs_ > 0, "n_eqs <= 0");
    static_assert(n_dims > 0, "n_dims <= 0");

    typedef real_t_ real_t; 
    typedef blitz::Array<real_t, n_dims> arr_t; 

    static const int n_eqs = n_eqs_;

    std::vector<int> n;
    arrvec_t<arr_t> C, psi[n_eqs];
    std::unordered_map<std::string, boost::ptr_vector<arrvec_t<arr_t>>> tmp; // intended for addressing with string(__FILE__)
    int span[n_dims];

    // accessor method for the Courant number field
    arr_t courant(int d = 0)  
    {   
      return C[d]; // TODO: what about halo? (in y)
    } 

    // ctor
    sharedmem_common()
      : n(n_eqs, 0)
    {}
  };

  template <int n_eqs = 1, typename real_t = float>
  struct sharedmem_1d : sharedmem_common<1, n_eqs, real_t>
  {
    blitz::Array<real_t, 1> state(int e)
    {
      return this->psi[e][ this->n[e] ](
	rng_t(0, this->span[0]-1)
      ).reindex({0});
    }

    // ctor
    sharedmem_1d(int s0)
    {
      this->span[0] = s0;
    }
  };

  template <int n_eqs = 1, typename real_t = float>
  struct sharedmem_2d : sharedmem_common<2, n_eqs, real_t>
  {
    blitz::Array<real_t, 2> state(int e)
    {
      return this->psi[e][ this->n[e] ](idx_t<2>({
	rng_t(0, this->span[0]-1),
	rng_t(0, this->span[1]-1)
      })).reindex({0, 0});
    }

    // ctor
    sharedmem_2d(int s0, int s1)
    {
      this->span[0] = s0;
      this->span[1] = s1;
    }
  };

  template <int n_eqs = 1, typename real_t = float>
  struct sharedmem_3d : sharedmem_common<3, n_eqs, real_t>
  {
    blitz::Array<real_t, 3> state(int e)
    {
      return this->psi[e][ this->n[e] ](idx_t<3>({
	rng_t(0, this->span[0]-1),
	rng_t(0, this->span[1]-1),
	rng_t(0, this->span[2]-1)
      })).reindex({0, 0, 0});
    }

    // ctor
    sharedmem_3d(int s0, int s1, int s2)
    {
      this->span[0] = s0;
      this->span[1] = s1;
      this->span[2] = s2;
    }
  };
}; // namespace advoocat
