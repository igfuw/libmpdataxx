/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <unordered_map>
#include <boost/ptr_container/ptr_vector.hpp>

#include "../../blitz.hpp"

namespace advoocat
{
  namespace concurr
  {
    namespace detail
    {
      template <
	typename real_t,
	int n_dims,
	int n_eqs
      >  
      class sharedmem_common
      {
	static_assert(n_eqs > 0, "n_eqs <= 0");
	static_assert(n_dims > 0, "n_dims <= 0");

        protected:

	int span[n_dims];

	public:

	std::vector<int> n;

	arrvec_t<blitz::Array<real_t, n_dims>> C, psi[n_eqs];

	std::unordered_map<
	  std::string, // intended for addressing with string(__FILE__)
	  boost::ptr_vector<arrvec_t<blitz::Array<real_t, n_dims>>>
	> tmp; 

	// accessor method for the Courant number field
	blitz::Array<real_t, n_dims> courant(int d = 0)  
	{   
	  return C[d]; // TODO: what about halo? (in y)
	}   

	// ctor
	sharedmem_common()
	  : n(n_eqs, 0)
	{}

        virtual void barrier()
        {
          assert(false && "sharedmem_common::barrier() called!");
        }
      };

      template<typename real_t, int n_dims, int n_eqs>
      class sharedmem
      {};

      template<typename real_t, int n_eqs>
      class sharedmem<real_t, 1, n_eqs> : public sharedmem_common<real_t, 1, n_eqs>
      {
	public:

	blitz::Array<real_t, 1> state(int e)
	{
	  return this->psi[e][ this->n[e] ](
	    rng_t(0, this->span[0]-1)
	  ).reindex({0});
	}

	// ctor
	sharedmem(int s0)
	{
	  this->span[0] = s0; 
	}
      };

      template<typename real_t, int n_eqs>
      class sharedmem<real_t, 2, n_eqs> : public sharedmem_common<real_t, 2, n_eqs>
      {
	public:

	blitz::Array<real_t, 2> state(int e)
	{
	  return this->psi[e][ this->n[e] ](idx_t<2>({
	    rng_t(0, this->span[0]-1),
	    rng_t(0, this->span[1]-1)
	  })).reindex({0, 0});
	}

	// ctor
	sharedmem(int s0, int s1)
	{
	  this->span[0] = s0; 
	  this->span[1] = s1; 
	}
      };

      template<typename real_t, int n_eqs>
      class sharedmem<real_t, 3, n_eqs> : public sharedmem_common<real_t, 3, n_eqs>
      {
	public:

	blitz::Array<real_t, 3> state(int e)
	{
	  return this->psi[e][ this->n[e] ](idx_t<3>({
	    rng_t(0, this->span[0]-1),
	    rng_t(0, this->span[1]-1),
	    rng_t(0, this->span[2]-1)
	  })).reindex({0, 0, 0});
	}

	// ctor
	sharedmem(int s0, int s1, int s2)
	{
	  this->span[0] = s0; 
	  this->span[1] = s1; 
	  this->span[2] = s2; 
	}
      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace advoocat
