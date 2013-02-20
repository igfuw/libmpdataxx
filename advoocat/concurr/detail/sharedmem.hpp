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
	int n_eqs,
        int n_tlev
      >  
      class sharedmem_common
      {
        using arr_t = blitz::Array<real_t, n_dims>; // TODO: use it in definitions of C, psi, tmp...

	static_assert(n_eqs > 0, "n_eqs <= 0");
	static_assert(n_dims > 0, "n_dims <= 0");
	static_assert(n_tlev > 0, "n_tlev <= 0");

        std::unique_ptr<blitz::Array<real_t, 1>> reducetmp;

        protected:

	int n, span[n_dims];

	public:

	arrvec_t<arr_t> C, psi[n_eqs];

	std::unordered_map< // TODO:! can string be used here (thread-safety!!!)
	  std::string, // intended for addressing with string(__FILE__)
	  boost::ptr_vector<arrvec_t<arr_t>>
	> tmp; 

	// accessor method for the Courant number field
	arr_t courant(int d = 0)  
	{   
	  return C[d]; // TODO: what about halo? (in y)
	}   

        virtual void barrier()
        {
          assert(false && "sharedmem_common::barrier() called!");
        }

        virtual int rank()
        {
          assert(false && "sharedmem_common::rank() called!");
        }
     
        void cycle()
        {
          if (rank() == 0) n = (n + 1) % n_tlev - n_tlev; // TODO: czy potrzebne - n_tlev?
        }

        // ctor
        sharedmem_common(int size)
        {
          reducetmp.reset(new blitz::Array<real_t, 1>(size));
          // TODO: fill reducetmp with NaNs
        }
  
        // concurrency-aware reductions
        real_t sum(const arr_t &arr)
        {
std::cerr << "summing ... " << std::endl;
          (*reducetmp)(rank()) = blitz::sum(arr); 
          barrier();
          real_t result = blitz::sum(*reducetmp);
std::ostringstream s;
s << "result: " << result << " rank=" << rank();
std::cerr << s << std::endl;
          barrier();
          return result;
        }

        real_t min(const arr_t &arr)
        {
          (*reducetmp)(rank()) = blitz::min(arr); 
          barrier();
          real_t result = blitz::min(*reducetmp);
          barrier();
          return result;
        }

        real_t max(const arr_t &arr)
        {
          (*reducetmp)(rank()) = blitz::max(arr); 
          barrier();
          real_t result = blitz::max(*reducetmp);
          barrier();
          return result;
        }
      };

      template<typename real_t, int n_dims, int n_eqs, int n_tlev>
      class sharedmem
      {};

      template<typename real_t, int n_eqs, int n_tlev>
      class sharedmem<real_t, 1, n_eqs, n_tlev> : public sharedmem_common<real_t, 1, n_eqs, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 1, n_eqs, n_tlev>;
	public:

	blitz::Array<real_t, 1> state(int e)
	{
	  return this->psi[e][ this->n ](
	    rng_t(0, this->span[0]-1)
	  ).reindex({0});
	}

	// ctor
	sharedmem(int s0, int size)
          : parent_t(size)
	{
	  this->span[0] = s0; 
	}
      };

      template<typename real_t, int n_eqs, int n_tlev>
      class sharedmem<real_t, 2, n_eqs, n_tlev> : public sharedmem_common<real_t, 2, n_eqs, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 2, n_eqs, n_tlev>;
	public:

	blitz::Array<real_t, 2> state(int e)
	{
	  return this->psi[e][ this->n ](idx_t<2>({
	    rng_t(0, this->span[0]-1),
	    rng_t(0, this->span[1]-1)
	  })).reindex({0, 0});
	}

	// ctor
	sharedmem(int s0, int s1, int size)
          : parent_t(size)
	{
	  this->span[0] = s0; 
	  this->span[1] = s1; 
	}
      };

      template<typename real_t, int n_eqs, int n_tlev>
      class sharedmem<real_t, 3, n_eqs, n_tlev> : public sharedmem_common<real_t, 3, n_eqs, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 3, n_eqs, n_tlev>;
	public:

	blitz::Array<real_t, 3> state(int e)
	{
	  return this->psi[e][ this->n ](idx_t<3>({
	    rng_t(0, this->span[0]-1),
	    rng_t(0, this->span[1]-1),
	    rng_t(0, this->span[2]-1)
	  })).reindex({0, 0, 0});
	}

	// ctor
	sharedmem(int s0, int s1, int s2, int size)
          : parent_t(size)
	{
	  this->span[0] = s0; 
	  this->span[1] = s1; 
	  this->span[2] = s2; 
	}
      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace advoocat
