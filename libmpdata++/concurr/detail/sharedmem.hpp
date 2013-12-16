/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <unordered_map>
#include <boost/ptr_container/ptr_vector.hpp>

#include <libmpdata++/blitz.hpp>

#include <array>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      template <
	typename real_t,
	int n_dims,
        int n_tlev
      >  
      class sharedmem_common
      {
        using arr_t = blitz::Array<real_t, n_dims>; 

	static_assert(n_dims > 0, "n_dims <= 0");
	static_assert(n_tlev > 0, "n_tlev <= 0");

        std::unique_ptr<blitz::Array<real_t, 1>> xtmtmp; // TODO: T_sumtype
        std::unique_ptr<blitz::Array<real_t, 2>> sumtmp; // TODO: T_sumtype

        protected:

        // for use in asserts only
        const int n_eqs;

	public:

	int n = 0;
        std::array<int, n_dims> span; 
        bool panic = false; // for multi-threaded SIGTERM handling

        // TODO: these are public because used from outside in alloc - could friendship help?
	arrvec_t<arr_t> GC;
        std::vector<arrvec_t<arr_t>> psi;
	std::unique_ptr<arr_t> G;

	std::unordered_map< 
	  const char*, // intended for addressing with __FILE__
	  boost::ptr_vector<arrvec_t<arr_t>>
	> tmp; 

        virtual void barrier()
        {
          assert(false && "sharedmem_common::barrier() called!");
        }

        virtual int rank()
        {
          assert(false && "sharedmem_common::rank() called!");
          throw;
        }
     
        void cycle()
        {
          barrier();
          if (rank() == 0) n = (n + 1) % n_tlev - n_tlev; // - n_tlev assures Python-type end-of-array cyclic behaviour works
          barrier();
        }

        // ctors
        // TODO: fill reducetmp with NaNs (or use 1-element arrvec_t - it's NaN-filled by default)
        sharedmem_common(const int &n_eqs, const std::array<int, 1> &span, const int &size)
          : n_eqs(n_eqs), n(0), psi(n_eqs), span(span) // TODO: is n(0) needed?
        {
          if (size > span[0]) throw std::exception(); // TODO: error_macro? / assert()?
          //sumtmp.reset(new blitz::Array<real_t, 2>(s0, 1));  // TODO: write a different sum that would't use sumtmp
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }

        sharedmem_common(const int &n_eqs, const std::array<int, 2> &span, const int &size)
          : n_eqs(n_eqs), n(0), psi(n_eqs), span(span)
        {
          if (size > span[0]) throw std::exception(); // TODO: error_macro?
          sumtmp.reset(new blitz::Array<real_t, 2>(span[0], 1));
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }

        sharedmem_common(const int &n_eqs, const std::array<int, 3> &span, const int &size)
          : n_eqs(n_eqs), n(0), psi(n_eqs), span(span)
        {
          if (size > span[0]) throw std::exception(); // TODO: error_macro?
          sumtmp.reset(new blitz::Array<real_t, 2>(span[0], span[1]));
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }
  
        /// @brief concurrency-aware 2D summation of array elements
        real_t sum(const arr_t &arr, const rng_t &i, const rng_t &j) // TODO: that's just for 2D
        {
	  // doing a two-step sum to reduce numerical error 
	  // and make parallel results reproducible
	  for (int c = i.first(); c <= i.last(); ++c) // TODO: optimise for i.count() == 1
          {
	    (*sumtmp)(c, 0) = blitz::kahan_sum(arr(c, j));
          }
          barrier();
          real_t result = blitz::kahan_sum(*sumtmp);
          barrier();
          return result;
        }

        /// @brief concurrency-aware 2D summation of a (element-wise) product of two arrays
        real_t sum(const arr_t &arr1, const arr_t &arr2, const rng_t &i, const rng_t &j) // TODO: that's just for 2D
        {
	  // doing a two-step sum to reduce numerical error 
	  // and make parallel results reproducible
	  for (int c = i.first(); c <= i.last(); ++c)
          {
	    (*sumtmp)(c, 0) = blitz::kahan_sum(arr1(c, j) * arr2(c, j)); 
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

        // this hack is introduced to allow to use neverDeleteData
        // and hence to not use BZ_THREADSAFE
        private:
        boost::ptr_vector<arr_t> tobefreed;
      
        public:
        arr_t *old(arr_t *arg)
        {
          tobefreed.push_back(arg);
          arr_t *ret = new arr_t(arg->dataFirst(), arg->shape(), blitz::neverDeleteData);
          ret->reindexSelf(arg->base());
          return ret;
        }
      };

      template<typename real_t, int n_dims, int n_tlev>
      class sharedmem
      {};

      template<typename real_t, int n_tlev>
      class sharedmem<real_t, 1, n_tlev> : public sharedmem_common<real_t, 1, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 1, n_tlev>;
        using parent_t::parent_t; // inheriting ctors

	public:

	// accessor methods
	blitz::Array<real_t, 1> advectee(int e = 0)
	{
          assert(e < this->n_eqs);
          assert(this->n < n_tlev);

          // returning just the domain interior, i.e. without halos
          // reindexing so that element 0 is at 0
	  return this->psi[e][ this->n ](
	    rng_t(0, this->span[0]-1)
	  ).reindex({0});
	}

	blitz::Array<real_t, 1> advector(int d = 0)  
	{   
          assert(d == 0);
          // returning the whole array but with the dimensiones
          // reindexed to make it more intuitive when working with index placeholders
          // (i.e. border between cell -1 and cell 0 is indexed with 0)
	  return this->GC[d](
            rng_t::all()
          ).reindex({this->GC[d].base(0)+1});
	}   

        blitz::Array<real_t, 1> g_factor()
        {
          // the same logic as in advectee() - see above
          return (*this->G)(
            rng_t(0, this->span[0]-1)
          ).reindex({0});
        }

      };

      template<typename real_t, int n_tlev>
      class sharedmem<real_t, 2, n_tlev> : public sharedmem_common<real_t, 2, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 2, n_tlev>;
        using parent_t::parent_t; // inheriting ctors

	public:

	blitz::Array<real_t, 2> advectee(int e = 0)
	{
          assert(e < this->n_eqs);
          assert(this->n < n_tlev);

	  return this->psi[e][ this->n ](idx_t<2>({
	    rng_t(0, this->span[0]-1),
	    rng_t(0, this->span[1]-1)
	  })).reindex({0, 0});
	}

	blitz::Array<real_t, 2> advector(int d = 0)  
	{   
          assert(d == 0 || d== 1);
          // returning the whole array (i.e. incl. haloes) but with the dimensiones
          // reindexed to make it more intuitive when working with index placeholders
          const rng_t all = rng_t::all();
          switch (d)
          { 
            case 0: return this->GC[d](all, all).reindex({this->GC[d].base(0)+1,this->GC[d].base(1)  }); 
            case 1: return this->GC[d](all, all).reindex({this->GC[d].base(0),  this->GC[d].base(1)+1}); 
            default: assert(false); throw;
          }
	}   

        blitz::Array<real_t, 2> g_factor()
        {
          // the same logic as in advectee() - see above
          return (*this->G)(idx_t<2>({
            rng_t(0, this->span[0]-1),
            rng_t(0, this->span[1]-1),
          })).reindex({0, 0});
        }

      };

      template<typename real_t, int n_tlev>
      class sharedmem<real_t, 3, n_tlev> : public sharedmem_common<real_t, 3, n_tlev>
      {
        using parent_t = sharedmem_common<real_t, 3, n_tlev>;
        using parent_t::parent_t; // inheriting ctors

	public:

	blitz::Array<real_t, 3> advectee(int e = 0)
	{
          assert(e < this->n_eqs);
          assert(this->n < n_tlev);

	  return this->psi[e][ this->n ](idx_t<3>({
	    rng_t(0, this->span[0]-1),
	    rng_t(0, this->span[1]-1),
	    rng_t(0, this->span[2]-1)
	  })).reindex({0, 0, 0});
	}

	blitz::Array<real_t, 3> advector(int d = 0)  
	{   
          assert(d == 0 || d == 1 || d == 2);
          // returning the whole array but with the dimensiones
          // reindexed to make it more intuitive when working with index placeholders
          const rng_t all = rng_t::all();
          switch (d)
          { 
            case 0: return this->GC[d](all, all, all).reindex({this->GC[d].base(0)+1,this->GC[d].base(1),  this->GC[d].base(2)  });  
            case 1: return this->GC[d](all, all, all).reindex({this->GC[d].base(0),  this->GC[d].base(1)+1,this->GC[d].base(2)  }); 
            case 2: return this->GC[d](all, all, all).reindex({this->GC[d].base(0),  this->GC[d].base(1),  this->GC[d].base(2)+1}); 
            default: assert(false); throw;
          }
	}   

        blitz::Array<real_t, 3> g_factor()
        {
          // the same logic as in advectee() - see above
          return (*this->G)(idx_t<3>({
            rng_t(0, this->span[0]-1),
            rng_t(0, this->span[1]-1),
            rng_t(0, this->span[2]-1)
          })).reindex({0, 0, 0});
        }

      };
    }; // namespace detail
  }; // namespace concurr
}; // namespace libmpdataxx
