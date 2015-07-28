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
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t>;
      using ix = typename ct_params_t::ix;

      protected:

      // member fields
      const rng_t im;

      void fill_stash() final
      {
        this->fill_stash_helper(0, ix::vip_i);
      }

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;

	if (!this->mem->G)
	{
	  this->mem->GC[0](im + h) = this->dt / this->di * .5 * (
	    this->stash[0](im    ) + 
	    this->stash[0](im + 1)
	  );
	} 
	else
	{ 
	  assert(false); // TODO: and if G is not const...
	}
      }

      void extrapolate_in_time() final
      {
	this->extrp(0, ix::vip_i);     
	this->xchng_sclr(this->stash[0]);      // filling halos 
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
      } 
    };

    // 2D version
    template <class ct_params_t> 
    class mpdata_rhs_vip<
      ct_params_t, 
      typename std::enable_if<ct_params_t::n_dims == 2>::type
    > : public detail::mpdata_rhs_vip_common<ct_params_t>
    {
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t>;
      using ix = typename ct_params_t::ix;

      protected:

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
	  this->mem->GC[d](pi<d>(i+h,j)) = this->dt / this->di * .5 * (
	    psi(pi<d>(i,    j)) + 
	    psi(pi<d>(i + 1,j))
	  );
	} 
	else
	{ 
	  assert(false); // TODO (perhaps better change definition of stash?)
	}
      }  

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;

	intrp<0>(this->stash[0], im, this->j^this->halo, this->di);
	intrp<1>(this->stash[1], jm, this->i^this->halo, this->dj);
        this->xchng_vctr_alng(this->mem->GC);
      }

      void extrapolate_in_time() final
      {
        using namespace libmpdataxx::arakawa_c; 

	this->extrp(0, ix::vip_i);     
	this->xchng_sclr(this->stash[0], this->i^this->halo, this->j^this->halo);      // filling halos 
	this->extrp(1, ix::vip_j);
	this->xchng_sclr(this->stash[1], this->i^this->halo, this->j^this->halo);      // filling halos 
      }
      
      void add_absorber() final
      {
	  typename ct_params_t::real_t factor = (ct_params_t::vip_vab == impl ? 0.5 : 1.);
          this->state(ix::vip_i)(this->ijk) -= factor * this->dt * 
            (*this->mem->vab_coeff)(this->ijk) * (this->state(ix::vip_i)(this->ijk) - this->mem->vab_relax[0](this->ijk));
          this->state(ix::vip_j)(this->ijk) -= factor * this->dt * 
            (*this->mem->vab_coeff)(this->ijk) * (this->state(ix::vip_j)(this->ijk) - this->mem->vab_relax[1](this->ijk));
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
      }
      
      static void alloc(typename parent_t::mem_t *mem, const typename parent_t::rt_params_t &p)
      {
	parent_t::alloc(mem, p);

        // allocate velocity absorber
        if (ct_params_t::vip_vab != 0)
        {
          mem->vab_coeff.reset(mem->old(new typename parent_t::arr_t(
                  parent_t::rng_sclr(p.grid_size[0]),
                  parent_t::rng_sclr(p.grid_size[1])
          )));
          
          for (int n = 0; n < ct_params_t::n_dims; ++n)
            mem->vab_relax.push_back(mem->old(new typename parent_t::arr_t(
                    parent_t::rng_sclr(p.grid_size[0]),
                    parent_t::rng_sclr(p.grid_size[1])
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
      using parent_t = detail::mpdata_rhs_vip_common<ct_params_t>;
      using ix = typename ct_params_t::ix;

      protected:

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
	  this->mem->GC[d](pi<d>(i+h, j, k)) = this->dt / this->di * .5 * (
	    psi(pi<d>(i,     j, k)) + 
	    psi(pi<d>(i + 1, j, k))
	  );
	} 
	else
	{ 
	  assert(false); // TODO (perhaps better change definition of stash?)
	}
      }  

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;

	intrp<0>(this->stash[0], im, this->j^this->halo, this->k^this->halo, this->di);
	intrp<1>(this->stash[1], jm, this->k^this->halo, this->i^this->halo, this->dj);
	intrp<2>(this->stash[2], km, this->i^this->halo, this->j^this->halo, this->dk);
      }

      void extrapolate_in_time() final
      {
        using namespace libmpdataxx::arakawa_c; 

	this->extrp(0, ix::vip_i);     
	this->xchng_sclr(this->stash[0],
                         this->i^this->halo,
                         this->j^this->halo,
                         this->k^this->halo);      // filling halos 
	this->extrp(1, ix::vip_j);
	this->xchng_sclr(this->stash[1],
                         this->i^this->halo,
                         this->j^this->halo,
                         this->k^this->halo);      // filling halos 
	this->extrp(2, ix::vip_k);
	this->xchng_sclr(this->stash[2],
                         this->i^this->halo,
                         this->j^this->halo,
                         this->k^this->halo);      // filling halos 
      }

      public:

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
      } 
    }; 
  }; // namespace solvers
}; // namespace libmpdataxx
