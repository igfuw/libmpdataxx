/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {

    // to be specialised
    template <typename ct_params_t, int minhalo = 0, class enableif = void>
    class mpdata_rhs_vip 
    {};

    // 1D version
    template <class ct_params_t, int minhalo> 
    class mpdata_rhs_vip<
      ct_params_t, minhalo,
      typename std::enable_if<ct_params_t::n_dims == 1>::type
    > : public detail::mpdata_rhs_vip_common<ct_params_t, minhalo>
    {
      using ix = typename ct_params_t::ix;

      protected:

      using solver_family = mpdata_rhs_vip_family_tag;
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t, minhalo>;

      // member fields
      const rng_t im;

      void interpolate_in_space(arrvec_t<typename parent_t::arr_t> &dst,
                                const arrvec_t<typename parent_t::arr_t> &src) final
      {
        using namespace libmpdataxx::arakawa_c;
	if (!this->mem->G)
	{
	  dst[0](im + h) = this->dt / this->di * .5 * (
	    src[0](im    ) + 
	    src[0](im + 1)
	  );
	} 
	else
	{ 
	  assert(false); // TODO: and if G is not const...
	}
        this->xchng_vctr_alng(dst, /*ad*/ false, /*cyclic*/ true);
      }

      void extrapolate_in_time() final
      {
	this->extrp(0, ix::vip_i);     
	this->xchng_sclr(this->vip_stash(0)[0]);      // filling halos 
      }

      public:

      // ctor
      mpdata_rhs_vip(
	typename parent_t::ctor_args_t args,
	const typename parent_t::rt_params_t &p
      ) : 
	parent_t(args, p),
	im(args.i.first() - 1, args.i.last())
      {
        assert(this->di != 0);

        this->vip_ixs = {ix::vip_i};
      } 
    };

    // 2D version
    template <class ct_params_t, int minhalo> 
    class mpdata_rhs_vip<
      ct_params_t, minhalo,
      typename std::enable_if<ct_params_t::n_dims == 2>::type
    > : public detail::mpdata_rhs_vip_common<ct_params_t, minhalo>
    {
      using ix = typename ct_params_t::ix;

      protected:

      using solver_family = mpdata_rhs_vip_family_tag;
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t, minhalo>;

      // member fields
      const rng_t im, jm;

      template<int d, class arr_t> 
      void intrp(
	arr_t &dst,
	const arr_t &src,
	const rng_t &i, 
	const rng_t &j, 
	const typename ct_params_t::real_t &di 
      )
      {   
	using idxperm::pi;
	using namespace arakawa_c;
  
	if (!this->mem->G)
	{
	  dst(pi<d>(i+h,j)) = this->dt / di * .5 * (
	    src(pi<d>(i,    j)) + 
	    src(pi<d>(i + 1,j))
	  );
	} 
	else
	{ 
	  dst(pi<d>(i+h,j)) = this->dt / di * .5 * (
	    (*this->mem->G)(pi<d>(i,    j)) * src(pi<d>(i,    j)) + 
	    (*this->mem->G)(pi<d>(i + 1,j)) * src(pi<d>(i + 1,j))
	  );
	}
      }  

      void interpolate_in_space(arrvec_t<typename parent_t::arr_t> &dst,
                                const arrvec_t<typename parent_t::arr_t> &src) final
      {
        using namespace arakawa_c; // for rng_t operator^
        auto ex = this->halo - 1;
	intrp<0>(dst[0], src[0], im^ex, this->j^ex, this->di);
	intrp<1>(dst[1], src[1], jm^ex, this->i^ex, this->dj);
        this->xchng_vctr_alng(dst, /*ad*/ false, /*cyclic*/ true);
        this->xchng_vctr_nrml(dst, this->ijk, ex, /*cyclic*/ false);
      }

      void extrapolate_in_time() final
      {
        using namespace libmpdataxx::arakawa_c; 

	this->extrp(0, ix::vip_i);
        // using xchng_pres because bcs have to be consistent with those used in
        // pressure solver to obtain non-divergent advector field
        auto ex = this->halo - 1;
	this->extrp(1, ix::vip_j);
	this->xchng_pres(this->vip_stash(0)[0], this->ijk, ex);
	this->xchng_pres(this->vip_stash(0)[1], this->ijk, ex);
      }

      public:
      
      // ctor
      mpdata_rhs_vip(
	typename parent_t::ctor_args_t args,
	const typename parent_t::rt_params_t &p
      ) : 
	parent_t(args, p),
	im(args.i.first() - 1, args.i.last()),
	jm(args.j.first() - 1, args.j.last())
      {
        assert(this->di != 0);
        assert(this->dj != 0);
        
        this->vip_ixs = {ix::vip_i, ix::vip_j};
      }
      
      static void alloc(
        typename parent_t::mem_t *mem, 
        const int &n_iters
      ) {
	parent_t::alloc(mem, n_iters);

        // allocate velocity absorber
        if (static_cast<vip_vab_t>(ct_params_t::vip_vab) != 0)
        {
          mem->vab_coeff.reset(mem->old(new typename parent_t::arr_t(
                  parent_t::rng_sclr(mem->grid_size[0]),
                  parent_t::rng_sclr(mem->grid_size[1])
          )));
          
          for (int n = 0; n < ct_params_t::n_dims; ++n)
            mem->vab_relax.push_back(mem->old(new typename parent_t::arr_t(
                    parent_t::rng_sclr(mem->grid_size[0]),
                    parent_t::rng_sclr(mem->grid_size[1])
            )));
        }
      }
    };
    
    // 3D version
    template <class ct_params_t, int minhalo> 
    class mpdata_rhs_vip<
      ct_params_t, minhalo,
      typename std::enable_if<ct_params_t::n_dims == 3>::type
    > : public detail::mpdata_rhs_vip_common<ct_params_t, minhalo>
    {
      using ix = typename ct_params_t::ix;

      protected:
      
      using solver_family = mpdata_rhs_vip_family_tag;
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t, minhalo>;

      // member fields
      const rng_t im, jm, km;

      template<int d, class arr_t> 
      void intrp(
	arr_t &dst,
	const arr_t &src,
	const rng_t &i, 
	const rng_t &j, 
	const rng_t &k, 
	const typename ct_params_t::real_t &di 
      )
      {   
	using idxperm::pi;
	using namespace arakawa_c;
  
	if (!this->mem->G)
	{
	  dst(pi<d>(i+h, j, k)) = this->dt / di * .5 * (
	    src(pi<d>(i,     j, k)) + 
	    src(pi<d>(i + 1, j, k))
	  );
	} 
	else
	{ 
	  dst(pi<d>(i+h, j, k)) = this->dt / di * .5 * (
	    (*this->mem->G)(pi<d>(i  , j, k)) * src(pi<d>(i,   j, k)) + 
	    (*this->mem->G)(pi<d>(i+1, j, k)) * src(pi<d>(i+1, j, k))
	  );
	}
      }  

      void interpolate_in_space(arrvec_t<typename parent_t::arr_t> &dst,
                                const arrvec_t<typename parent_t::arr_t> &src) final
      {
        using namespace arakawa_c; // for rng_t operator^
        auto ex = this->halo - 1;
	intrp<0>(dst[0], src[0], im^ex, this->j^ex, this->k^ex, this->di);
	intrp<1>(dst[1], src[1], jm^ex, this->k^ex, this->i^ex, this->dj);
	intrp<2>(dst[2], src[2], km^ex, this->i^ex, this->j^ex, this->dk);
        this->xchng_vctr_alng(dst, /*ad*/ false, /*cyclic*/ true);
        this->xchng_vctr_nrml(dst, this->ijk, ex, /*cyclic*/ false);
      }

      void extrapolate_in_time() final
      {
        using namespace libmpdataxx::arakawa_c; 

        // using xchng_pres because bcs have to be consistent with those used in
        // pressure solver to obtain non-divergent advector field
        auto ex = this->halo - 1;
	this->extrp(0, ix::vip_i);     
	this->extrp(1, ix::vip_j);
	this->extrp(2, ix::vip_k);
	this->xchng_pres(this->vip_stash(0)[0], this->ijk, ex);
	this->xchng_pres(this->vip_stash(0)[1], this->ijk, ex);
	this->xchng_pres(this->vip_stash(0)[2], this->ijk, ex);
      }

      public:
      
      static void alloc(
        typename parent_t::mem_t *mem, 
        const int &n_iters
      ) 
      {
	parent_t::alloc(mem, n_iters);

        // allocate velocity absorber
        if (static_cast<vip_vab_t>(ct_params_t::vip_vab) != 0)
        {
          mem->vab_coeff.reset(mem->old(new typename parent_t::arr_t(
                  parent_t::rng_sclr(mem->grid_size[0]),
                  parent_t::rng_sclr(mem->grid_size[1]),
                  parent_t::rng_sclr(mem->grid_size[2])
          )));
          
          for (int n = 0; n < ct_params_t::n_dims; ++n)
            mem->vab_relax.push_back(mem->old(new typename parent_t::arr_t(
                    parent_t::rng_sclr(mem->grid_size[0]),
                    parent_t::rng_sclr(mem->grid_size[1]),
                    parent_t::rng_sclr(mem->grid_size[2])
            )));
        }
      }

      // ctor
      mpdata_rhs_vip(
	typename parent_t::ctor_args_t args,
	const typename parent_t::rt_params_t &p
      ) : 
	parent_t(args, p),
	im(args.i.first() - 1, args.i.last()),
	jm(args.j.first() - 1, args.j.last()),
	km(args.k.first() - 1, args.k.last())
      {
        assert(this->di != 0);
        assert(this->dj != 0);
        assert(this->dk != 0);
        
        this->vip_ixs = {ix::vip_i, ix::vip_j, ix::vip_k};
      } 
    }; 
  } // namespace solvers
} // namespace libmpdataxx
