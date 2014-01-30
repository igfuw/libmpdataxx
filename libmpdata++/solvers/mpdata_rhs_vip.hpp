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
      const typename ct_params_t::real_t di;

      void fill_stash() final
      {
	if (ix::vip_den == -1) 
	  this->stash[0](this->ijk) = this->psi_n(ix::vip_i)(this->ijk);
	else if (this->initial_h_non_zero)
          this->stash[0](this->ijk) = this->psi_n(ix::vip_i)(this->ijk) / this->psi_n(ix::vip_den)(this->ijk);
        else
	{
	  this->stash[0](this->ijk) = where(
	    // if
	    this->psi_n(ix::vip_den)(this->ijk) > 0,
	    // then
	    this->psi_n(ix::vip_i)(this->ijk) / this->psi_n(ix::vip_den)(this->ijk),
	    // else
            0
	  );
	}
      }

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;

	if (!this->mem->G)
	{
	  this->mem->GC[0](im+h) = this->dt / di * .5 * (
	    this->stash[0](im    ) + 
	    this->stash[0](im + 1)
	  );
	} 
	else
	{ 
	  assert(false); // TODO
	}
      }

      void extrapolate_in_time() final
      {
        using namespace arakawa_c;

	this->stash[0](this->ijk) /= -2.;

	if (ix::vip_den == -1) 
	  this->stash[0](this->ijk) += 3./2 * this->psi_n(ix::vip_i)(this->ijk);
	else if (this->initial_h_non_zero)
	  this->stash[0](this->ijk) += 3./2 * this->psi_n(ix::vip_i)(this->ijk) / this->psi_n(ix::vip_den)(this->ijk);
	else
        {
          // TODO: hint_nozeroh
	  this->stash[0](this->ijk) += where(
            // if
            this->psi_n(ix::vip_den)(this->ijk) > 0,
            // then
            3./2 * this->psi_n(ix::vip_i)(this->ijk) / this->psi_n(ix::vip_den)(this->ijk),
            // else
            0
          ); 
        }

	this->xchng(this->stash[0]);      // filling halos 
      }

      public:

      struct rt_params_t : parent_t::rt_params_t 
      { 
	typename ct_params_t::real_t di=0;
      };

      // ctor
      mpdata_rhs_vip(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) : 
	parent_t(args, p),
	im(args.i.first() - 1, args.i.last()),
	di(p.di)
      {
        assert(di != 0);
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
      const typename ct_params_t::real_t di, dj;

      void fill_stash() final
      {
	if (ix::vip_den == -1) 
	{
	  this->stash[0](this->ijk) = this->psi_n(ix::vip_i)(this->ijk);
	  this->stash[1](this->ijk) = this->psi_n(ix::vip_j)(this->ijk);
	}
	else
	{
	  this->stash[0](this->ijk) = this->psi_n(ix::vip_i)(this->ijk) / this->psi_n(ix::vip_den)(this->ijk); // TODO: what if density == 0?
	  this->stash[1](this->ijk) = this->psi_n(ix::vip_j)(this->ijk) / this->psi_n(ix::vip_den)(this->ijk); // TODO: what if density == 0?
	} 
      }

      template <int d>
      void extrp(int e) // extrapolate velocity field in time to t+1/2
      {                 // (write the result to stash since we don't need previous state any more)
        using namespace arakawa_c;

	this->stash[d](this->ijk) /= -2.;

	if (ix::vip_den == -1) 
	  this->stash[d](this->ijk) += 3./2 * this->psi_n(e)(this->ijk);
	else
	  this->stash[d](this->ijk) += 3./2 * (this->psi_n(e)(this->ijk) / this->psi_n(ix::vip_den)(this->ijk)); // TODO: what if density == 0?

	this->xchng(this->stash[d], this->i^this->halo, this->j^this->halo);      // filling halos 
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
	  assert(false); // TODO
	}
      }  

      void interpolate_in_space() final
      {
        using namespace libmpdataxx::arakawa_c;

	intrp<0>(this->stash[0], im, this->j^this->halo, di);
	intrp<1>(this->stash[1], jm, this->i^this->halo, dj);
      }

      void extrapolate_in_time() final
      {
	extrp<0>(ix::vip_i);     
	extrp<1>(ix::vip_j);
      }

      public:

      struct rt_params_t : parent_t::rt_params_t 
      { 
	typename ct_params_t::real_t di=0, dj=0;
      };

      // ctor
      mpdata_rhs_vip(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) : 
	parent_t(args, p),
	im(args.i.first() - 1, args.i.last()),
	jm(args.j.first() - 1, args.j.last()),
	di(p.di), dj(p.dj)
      {
        assert(di != 0);
        assert(dj != 0);
      } 
    }; 
  }; // namespace solvers
}; // namespace libmpdataxx
