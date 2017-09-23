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
    template <typename ct_params_t, class enableif = void>
    class mpdata_rhs_vip 
    {};

    // 1D version
    template <class ct_params_t> 
    class mpdata_rhs_vip<
      ct_params_t,
      typename std::enable_if<ct_params_t::n_dims == 1>::type
    > : public detail::mpdata_rhs_vip_common<ct_params_t>
    {
      using ix = typename ct_params_t::ix;

      protected:

      using solver_family = mpdata_rhs_vip_family_tag;
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t>;

      // member fields
      const rng_t im;

      void fill_stash() final
      {
        this->fill_stash_helper(0, ix::vip_i);
      }

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;
        
        const auto off = ct_params_t::var_dt ? ct_params_t::n_dims : 0;

	if (!this->mem->G)
	{
	  this->mem->GC[0](im + h) = this->dt / this->di * .5 * (
	    this->stash[0 + off](im    ) + 
	    this->stash[0 + off](im + 1)
	  );
	} 
	else
	{ 
	  assert(false); // TODO: and if G is not const...
	}
      }

      void extrapolate_in_time() final
      {
        const auto off = ct_params_t::var_dt ? ct_params_t::n_dims : 0;
	this->extrp(0, ix::vip_i);     
	this->xchng_sclr(this->stash[0 + off]);      // filling halos 
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
    template <class ct_params_t> 
    class mpdata_rhs_vip<
      ct_params_t, 
      typename std::enable_if<ct_params_t::n_dims == 2>::type
    > : public detail::mpdata_rhs_vip_common<ct_params_t>
    {
      using ix = typename ct_params_t::ix;

      protected:

      using solver_family = mpdata_rhs_vip_family_tag;
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t>;

      // member fields
      const rng_t im, jm;

      void fill_stash() final 
      {
        this->fill_stash_helper(0, ix::vip_i);
        this->fill_stash_helper(1, ix::vip_j);
      }

      template<int d, class arr_t> 
      void intrp(
	const arr_t psi,
	const rng_t &i, 
	const rng_t &j, 
	const typename ct_params_t::real_t &di 
      )
      {   
	using idxperm::pi;
	using namespace arakawa_c;
  
	if (!this->mem->G)
	{
	  this->mem->GC[d](pi<d>(i+h,j)) = this->dt / di * .5 * (
	    psi(pi<d>(i,    j)) + 
	    psi(pi<d>(i + 1,j))
	  );
	} 
	else
	{ 
	  this->mem->GC[d](pi<d>(i+h,j)) = this->dt / di * .5 * (
	    (*this->mem->G)(pi<d>(i,    j)) * psi(pi<d>(i,    j)) + 
	    (*this->mem->G)(pi<d>(i + 1,j)) * psi(pi<d>(i + 1,j))
	  );
	}
      }  

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;

        const auto off = ct_params_t::var_dt ? ct_params_t::n_dims : 0;

	intrp<0>(this->stash[0 + off], im, this->j^this->halo, this->di);
	intrp<1>(this->stash[1 + off], jm, this->i^this->halo, this->dj);
        this->xchng_vctr_alng(this->mem->GC);
        auto ex = this->halo - 1;
        this->xchng_vctr_nrml(this->mem->GC, this->i^ex, this->j^ex);
      }

      void extrapolate_in_time() final
      {
        using namespace libmpdataxx::arakawa_c; 

        const auto off = ct_params_t::var_dt ? ct_params_t::n_dims : 0;

	this->extrp(0, ix::vip_i);     
	this->xchng_sclr(this->stash[0 + off], this->i^this->halo, this->j^this->halo);      // filling halos 
	this->extrp(1, ix::vip_j);
	this->xchng_sclr(this->stash[1 + off], this->i^this->halo, this->j^this->halo);      // filling halos 
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
    template <class ct_params_t> 
    class mpdata_rhs_vip<
      ct_params_t, 
      typename std::enable_if<ct_params_t::n_dims == 3>::type
    > : public detail::mpdata_rhs_vip_common<ct_params_t>
    {
      using ix = typename ct_params_t::ix;

      protected:
      
      using solver_family = mpdata_rhs_vip_family_tag;
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t>;

      // member fields
      const rng_t im, jm, km;

      void fill_stash() final 
      {
        this->fill_stash_helper(0, ix::vip_i);
        this->fill_stash_helper(1, ix::vip_j);
        this->fill_stash_helper(2, ix::vip_k);
      }

      template<int d, class arr_t> 
      void intrp(
	const arr_t psi,
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
	  this->mem->GC[d](pi<d>(i+h, j, k)) = this->dt / di * .5 * (
	    psi(pi<d>(i,     j, k)) + 
	    psi(pi<d>(i + 1, j, k))
	  );
	} 
	else
	{ 
	  this->mem->GC[d](pi<d>(i+h, j, k)) = this->dt / di * .5 * (
	    (*this->mem->G)(pi<d>(i  , j, k)) * psi(pi<d>(i,   j, k)) + 
	    (*this->mem->G)(pi<d>(i+1, j, k)) * psi(pi<d>(i+1, j, k))
	  );
	}
      }  

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;
        
        const auto off = ct_params_t::var_dt ? ct_params_t::n_dims : 0;

	intrp<0>(this->stash[0 + off], im, this->j^this->halo, this->k^this->halo, this->di);
	intrp<1>(this->stash[1 + off], jm, this->k^this->halo, this->i^this->halo, this->dj);
	intrp<2>(this->stash[2 + off], km, this->i^this->halo, this->j^this->halo, this->dk);
        this->xchng_vctr_alng(this->mem->GC);
        auto ex = this->halo - 1;
        this->xchng_vctr_nrml(this->mem->GC, this->i^ex, this->j^ex, this->k^ex);
      }

      void extrapolate_in_time() final
      {
        using namespace libmpdataxx::arakawa_c; 

        const auto off = ct_params_t::var_dt ? ct_params_t::n_dims : 0;

	this->extrp(0, ix::vip_i);     
	this->xchng_sclr(this->stash[0 + off],
                         this->i^this->halo,
                         this->j^this->halo,
                         this->k^this->halo);      // filling halos 
	this->extrp(1, ix::vip_j);
	this->xchng_sclr(this->stash[1 + off],
                         this->i^this->halo,
                         this->j^this->halo,
                         this->k^this->halo);      // filling halos 
	this->extrp(2, ix::vip_k);
	this->xchng_sclr(this->stash[2 + off],
                         this->i^this->halo,
                         this->j^this->halo,
                         this->k^this->halo);      // filling halos 
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
