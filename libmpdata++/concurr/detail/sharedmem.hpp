/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <unordered_map>
#include <boost/ptr_container/ptr_vector.hpp>

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/formulae/arakawa_c.hpp>

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

        // TODO: T_sumtype (perhaps worh using double even if summing floats?)
        std::unique_ptr<blitz::Array<real_t, 1>> xtmtmp; 
        std::unique_ptr<blitz::Array<real_t, 1>> sumtmp;

        protected:

        blitz::TinyVector<int, n_dims> origin;

	public:

	int n = 0;
	const int size;
        std::array<rng_t, n_dims> grid_size; 
        bool panic = false; // for multi-threaded SIGTERM handling

        // TODO: these are public because used from outside in alloc - could friendship help?
	arrvec_t<arr_t> GC, ndt_GC, ndtt_GC;
        std::vector<arrvec_t<arr_t>> psi; // TODO: since n_eqns is known, could make it an std::array!
	std::unique_ptr<arr_t> G;
	std::unique_ptr<arr_t> vab_coeff; // velocity absorber coefficient
	arrvec_t<arr_t> vab_relax; // velocity absorber relaxed state
        arrvec_t<arr_t> khn_tmp; // Kahan sum for donor-cell

	std::unordered_map< 
	  const char*, // intended for addressing with __FILE__
	  boost::ptr_vector<arrvec_t<arr_t>>
	> tmp;
        
        // list of temporary fields that can be accessed from outside of concurr
	std::unordered_map< 
          std::string,
          std::pair<const char*, int>
        > avail_tmp;

        virtual void barrier()
        {
          assert(false && "sharedmem_common::barrier() called!");
        }
     
        void cycle(const int &rank)
        {
          barrier();
          if (rank == 0) n = (n + 1) % n_tlev - n_tlev; // - n_tlev assures Python-type end-of-array cyclic behaviour works
          barrier();
        }

        // ctors
        // TODO: fill reducetmp with NaNs (or use 1-element arrvec_t - it's NaN-filled by default)
        sharedmem_common(const std::array<int, n_dims> &grid_size, const int &size)
          : n(0), size(size) // TODO: is n(0) needed?
        {
          for (int d = 0; d < n_dims; ++d) 
          {
            this->grid_size[d] = rng_t(0, grid_size[d]-1);
            origin[d] = this->grid_size[d].first();
          }

          if (size > grid_size[0]) 
            throw std::runtime_error("number of subdomains greater than number of gridpoints");

          if (n_dims != 1) 
            sumtmp.reset(new blitz::Array<real_t, 1>(grid_size[0]));
          xtmtmp.reset(new blitz::Array<real_t, 1>(size));
        }

        /// @brief concurrency-aware summation of array elements
        real_t sum(const arr_t &arr, const idx_t<n_dims> &ijk, const bool sum_khn)
        {
	  // doing a two-step sum to reduce numerical error 
	  // and make parallel results reproducible
	  for (int c = ijk[0].first(); c <= ijk[0].last(); ++c) // TODO: optimise for i.count() == 1
          {
            auto slice_idx = ijk;
            slice_idx.lbound(0) = c;
            slice_idx.ubound(0) = c;

            if (sum_khn)
	      (*sumtmp)(c) = blitz::kahan_sum(arr(slice_idx));
            else
	      (*sumtmp)(c) = blitz::sum(arr(slice_idx));
          }
          barrier();
          real_t result;
          if (sum_khn)
            result = blitz::kahan_sum(*sumtmp);
          else
            result = blitz::sum(*sumtmp);
          barrier();
          return result;
        }

        /// @brief concurrency-aware summation of a (element-wise) product of two arrays
        real_t sum(const arr_t &arr1, const arr_t &arr2, const idx_t<n_dims> &ijk, const bool sum_khn)
        {
	  // doing a two-step sum to reduce numerical error 
	  // and make parallel results reproducible
	  for (int c = ijk[0].first(); c <= ijk[0].last(); ++c)
          {
            auto slice_idx = ijk;
            slice_idx.lbound(0) = c;
            slice_idx.ubound(0) = c;

            if (sum_khn)
	      (*sumtmp)(c) = blitz::kahan_sum(arr1(slice_idx) * arr2(slice_idx));
            else
	      (*sumtmp)(c) = blitz::sum(arr1(slice_idx) * arr2(slice_idx)); 
          }
          barrier();
          real_t result;
          if (sum_khn)
            result = blitz::kahan_sum(*sumtmp);
          else
            result = blitz::sum(*sumtmp);
          barrier();
          return result;
        }

        real_t min(const int &rank, const arr_t &arr)
        {
          (*xtmtmp)(rank) = blitz::min(arr); 
          barrier();
          real_t result = blitz::min(*xtmtmp);
          barrier();
          return result;
        }

        real_t max(const int &rank, const arr_t &arr)
        {
          (*xtmtmp)(rank) = blitz::max(arr); 
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
        arr_t *never_delete(arr_t *arg)
        {
          arr_t *ret = new arr_t(arg->dataFirst(), arg->shape(), blitz::neverDeleteData);
          ret->reindexSelf(arg->base());
          return ret;
        }

        arr_t *old(arr_t *arg)
        {
          tobefreed.push_back(arg);
          arr_t *ret = never_delete(arg);
          return ret;
        }

        private:
        // helper methods to define subdomain ranges
        static int min(const int &span, const int &rank, const int &size) 
        {
          return rank * span / size; 
        }

        static int max(const int &span, const int &rank, const int &size) 
        {
          return min(span, rank + 1, size) - 1;  
        }

        public:
        static rng_t slab(
          const rng_t &span,
          const int &rank = 0,  
          const int &size = 1 
        ) {
          return rng_t(
            span.first() + min(span.length(), rank, size),
            span.first() + max(span.length(), rank, size)
          );
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
          assert(this->n < n_tlev);

          // returning just the domain interior, i.e. without halos
          // reindexing so that element 0 is at 0
	  return this->psi[e][ this->n ](
            this->grid_size[0]
	  ).reindex(this->origin);
	}

	blitz::Array<real_t, 1> advector(int d = 0)  
	{   
          using namespace arakawa_c;
          assert(d == 0);
          // returning just the domain interior, i.e. without halos
          // reindexed to make it more intuitive when working with index placeholders
          // (i.e. border between cell 0 and cell 1 is indexed with 0)
	  return this->GC[d](
            this->grid_size[0]^(-1)^h
          ).reindex(this->origin);
	}   

        blitz::Array<real_t, 1> g_factor()
        {
          // a sanity check
          if (this->G.get() == nullptr) 
            throw std::runtime_error("g_factor() called with nug option unset?");

          // the same logic as in advectee() - see above
          return (*this->G)(
            this->grid_size[0]
          ).reindex(this->origin);
        }
        
        blitz::Array<real_t, 1> vab_coefficient()
        {
          throw std::logic_error("absorber not yet implemented in 1d");
        }
	
        blitz::Array<real_t, 1> vab_relaxed_state(int d = 0)  
	{   
          throw std::logic_error("absorber not yet implemented in 1d");
	}   
        
        blitz::Array<real_t, 1> sclr_array(const std::string& name, int n = 0)
	{
          return this->tmp.at(this->avail_tmp[name].first)[this->avail_tmp[name].second][n](
            this->grid_size[0]
          ).reindex(this->origin);
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
          assert(this->n < n_tlev);

	  return this->psi[e][ this->n ](
	    this->grid_size[0],
	    this->grid_size[1]
	  ).reindex(this->origin);
	}

	blitz::Array<real_t, 2> advector(int d = 0)  
	{   
          using namespace arakawa_c;
          assert(d == 0 || d== 1);
          // returning just the domain interior, i.e. without halos
          // reindexed to make it more intuitive when working with index placeholders
          switch (d)
          { 
            case 0: return this->GC[d](this->grid_size[0]^(-1)^h, this->grid_size[1]).reindex(this->origin); 
            case 1: return this->GC[d](this->grid_size[0], this->grid_size[1]^(-1)^h).reindex(this->origin); 
            default: assert(false); throw;
          }
	}   

        blitz::Array<real_t, 2> g_factor()
        {
          // a sanity check
          if (this->G.get() == nullptr) 
            throw std::runtime_error("g_factor() called with nug option unset?");

          // the same logic as in advectee() - see above
          return (*this->G)(
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
        }
        
        blitz::Array<real_t, 2> vab_coefficient()
        {
          // a sanity check
          if (this->vab_coeff.get() == nullptr) 
            throw std::runtime_error("vab_coeff() called with option vip_vab unset?");

          // the same logic as in advectee() - see above
          return (*this->vab_coeff)(
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
        }
	
        blitz::Array<real_t, 2> vab_relaxed_state(int d = 0)  
	{   
          assert(d == 0 || d== 1);
          // a sanity check
          if (this->vab_coeff.get() == nullptr) 
            throw std::runtime_error("vab_relaxed_state() called with option vip_vab unset?");
          // the same logic as in advectee() - see above
	  return this->vab_relax[d](
	    this->grid_size[0],
	    this->grid_size[1]
	  ).reindex(this->origin);
	}   
        
        blitz::Array<real_t, 2> sclr_array(const std::string& name, int n = 0)
	{
          return this->tmp.at(this->avail_tmp[name].first)[this->avail_tmp[name].second][n](
            this->grid_size[0],
            this->grid_size[1]
          ).reindex(this->origin);
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
          assert(this->n < n_tlev);

	  return this->psi[e][ this->n ](
	    this->grid_size[0],
	    this->grid_size[1],
	    this->grid_size[2]
	  ).reindex(this->origin);
	}

	blitz::Array<real_t, 3> advector(int d = 0)  
	{   
          using namespace arakawa_c;
          assert(d == 0 || d == 1 || d == 2);
          // returning just the domain interior, i.e. without halos
          // reindexed to make it more intuitive when working with index placeholders
          switch (d)
          { 
            case 0: return this->GC[d](this->grid_size[0]^(-1)^h,
                                       this->grid_size[1],
                                       this->grid_size[2]).reindex(this->origin);  
            case 1: return this->GC[d](this->grid_size[0],
                                       this->grid_size[1]^(-1)^h,
                                       this->grid_size[2]).reindex(this->origin);  
            case 2: return this->GC[d](this->grid_size[0],
                                       this->grid_size[1],
                                       this->grid_size[2]^(-1)^h).reindex(this->origin);  
            default: assert(false); throw;
          }
	}   

        blitz::Array<real_t, 3> g_factor()
        {
          // a sanity check
          if (this->G.get() == nullptr) 
            throw std::runtime_error("g_factor() called with nug option unset?");

          // the same logic as in advectee() - see above
          return (*this->G)(
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
        }

        blitz::Array<real_t, 3> vab_coefficient()
        {
          // a sanity check
          if (this->vab_coeff.get() == nullptr) 
            throw std::runtime_error("vab_coeff() called with option vip_vab unset?");

          // the same logic as in advectee() - see above
          return (*this->vab_coeff)(
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
        }
	
        blitz::Array<real_t, 3> vab_relaxed_state(int d = 0)  
	{   
          assert(d == 0 || d == 1 || d == 2);
          // a sanity check
          if (this->vab_coeff.get() == nullptr) 
            throw std::runtime_error("vab_relaxed_state() called with option vip_vab unset?");
          // the same logic as in advectee() - see above
	  return this->vab_relax[d](
	    this->grid_size[0],
	    this->grid_size[1],
            this->grid_size[2]
	  ).reindex(this->origin);
	}   

        blitz::Array<real_t, 3> sclr_array(const std::string& name, int n = 0)
	{
          return this->tmp.at(this->avail_tmp[name].first)[this->avail_tmp[name].second][n](
            this->grid_size[0],
            this->grid_size[1],
            this->grid_size[2]
          ).reindex(this->origin);
	}
      };
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
