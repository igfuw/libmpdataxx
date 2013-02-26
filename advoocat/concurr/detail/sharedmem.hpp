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

        std::unique_ptr<blitz::Array<real_t, 1>> xtmtmp; // TODO: T_sumtype
        std::unique_ptr<blitz::Array<real_t, 2>> sumtmp; // TODO: T_sumtype

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
          barrier();
          if (rank() == 0) n = (n + 1) % n_tlev - n_tlev; // TODO: czy potrzebne - n_tlev?
          barrier();
        }

        // ctors
        // TODO: fill reducetmp with NaNs
        sharedmem_common(int s0, int size)
        {
          if (size > s0) throw std::exception(); // TODO: error_macro?
          //sumtmp.reset(new blitz::Array<real_t, 2>(s0, 1));  // TODO: write a different sum that would't use sumtmp
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }
        sharedmem_common(int s0, int s1, int size)
        {
          if (size > s0) throw std::exception(); // TODO: error_macro?
          sumtmp.reset(new blitz::Array<real_t, 2>(s0, 1));
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }
        sharedmem_common(int s0, int s1, int s2, int size)
        {
          if (size > s0) throw std::exception(); // TODO: error_macro?
          sumtmp.reset(new blitz::Array<real_t, 2>(s0, s1));
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }
  
        // concurrency-aware reductions
        real_t sum(const arr_t &arr, const rng_t &i, const rng_t &j)
        {
          {
            // doing a two-step sum to reduce numerical error 
            // and make parallel results reproducible
            for (int c = i.first(); c <= i.last(); ++c)
              (*sumtmp)(c, 0) = blitz::kahan_sum(arr(c, j));
          }
          barrier();
          real_t result = blitz::kahan_sum(*sumtmp);
          barrier();
          return result;
        }

        real_t min(const arr_t &arr)
        {
          (*xtmtmp)(rank()) = blitz::min(arr); 
          barrier();
          real_t result = blitz::min(*xtmtmp);
          barrier();
          return result;
        }

        real_t max(const arr_t &arr)
        {
          (*xtmtmp)(rank()) = blitz::max(arr); 
          barrier();
          real_t result = blitz::max(*xtmtmp);
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
          : parent_t(s0, size)
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
          : parent_t(s0, s1, size)
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
          : parent_t(s0, s1, s2, size)
	{
	  this->span[0] = s0; 
	  this->span[1] = s1; 
	  this->span[2] = s2; 
	}
      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace advoocat
