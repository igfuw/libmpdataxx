// 1D MPI ``remote'' boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/remote_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int dim>    
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     dim,
      typename std::enable_if<
        knd == remote && 
        dir == left   && 
        n_dims == 1
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {
      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = typename parent_t::arr_t;
      using idx_t = typename parent_t::idx_t;
      using idx_ctor_arg_t = blitz::TinyVector<rng_t, n_dims>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : -1;

      public:

      void fill_halos_sclr(arr_t &a, const bool deriv = false)
      {
        this->xchng(a, idx_t(idx_ctor_arg_t(this->left_intr_sclr + off)), idx_t(idx_ctor_arg_t(this->left_halo_sclr)));
      }

      void fill_halos_pres(arr_t &a)
      {
        fill_halos_sclr(a);
      }
      
      void save_edge_vel(const arr_t &) {}

      void set_edge_pres(arr_t &, int) {}

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const bool ad = false)
      {
        if(!this->is_cyclic) 
        {
          if(halo == 1)
            this->send(av[0], idx_t(idx_ctor_arg_t(this->left_intr_vctr + off))); 
          else
            // processes fill vectors to the left of their domain
            this->xchng(av[0], idx_t(idx_ctor_arg_t(this->left_intr_vctr + off)), idx_t(idx_ctor_arg_t((this->left_halo_vctr^h)^(-1)))); 
        }
        else 
          // cyclic should communicate both ways 
          this->xchng(av[0], idx_t(idx_ctor_arg_t(this->left_intr_vctr + off)), idx_t(idx_ctor_arg_t(this->left_halo_vctr))); 
      }
      
      void fill_halos_sgs_div(arr_t &a)
      {
        fill_halos_sclr(a);
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const int offset = 0)
      {
	using namespace idxperm;
        // the same logic as fill_halos_vctr_alng but have to consider offset ... TODO: find a way to reuse !
        if(!this->is_cyclic) 
        {
          if(halo == 1)
            this->send(av[0 + offset], idx_t(idx_ctor_arg_t(this->left_intr_vctr + off))); 
          else
            // processes fill vectors to the left of their domain
            this->xchng(av[0 + offset], idx_t(idx_ctor_arg_t(this->left_intr_vctr + off)), idx_t(idx_ctor_arg_t((this->left_halo_vctr^h)^(-1)))); 
        }
        else 
          // cyclic should communicate both ways 
          this->xchng(av[0 + offset], idx_t(idx_ctor_arg_t(this->left_intr_vctr + off)), idx_t(idx_ctor_arg_t(this->left_halo_vctr))); 
      }
      
      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const real_t)
      {
        fill_halos_vctr_alng(av);
      }

      void fill_halos_vctr_nrml(arr_t &a)                 
      {                                                                         
        fill_halos_sclr(a);                                                  
      }  
      
      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const bool ad = false)
      {
        fill_halos_vctr_alng(av, ad);
      }

      void fill_halos_vctr_nrml_cyclic(arr_t &a)
      {
        fill_halos_vctr_nrml(a);
      }
    };



    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int dim>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     dim,
      typename std::enable_if<
        knd == remote &&
        dir == rght   &&
        n_dims == 1
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {
      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = typename parent_t::arr_t;
      using idx_t = typename parent_t::idx_t;
      using idx_ctor_arg_t = blitz::TinyVector<rng_t, n_dims>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : 1;

      public:

      void fill_halos_sclr(arr_t &a, const bool deriv = false)
      {
        this->xchng(a, idx_t(idx_ctor_arg_t(this->rght_intr_sclr + off)), idx_t(idx_ctor_arg_t(this->rght_halo_sclr)));
      }

      void fill_halos_pres(arr_t &a)
      {
        fill_halos_sclr(a);
      }
      
      void save_edge_vel(const arr_t &) {}

      void set_edge_pres(arr_t &, int) {}


      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const bool ad = false)
      {
        if(!this->is_cyclic) 
        {
          if(halo == 1)
            this->recv(av[0],  idx_t(idx_ctor_arg_t(this->rght_halo_vctr))); 
          else
            this->xchng(av[0], idx_t(idx_ctor_arg_t(((this->rght_intr_vctr + off)^h)^(-1))), idx_t(idx_ctor_arg_t(this->rght_halo_vctr))); 
        }
        else
          this->xchng(av[0], idx_t(idx_ctor_arg_t(this->rght_intr_vctr + off)), idx_t(idx_ctor_arg_t(this->rght_halo_vctr))); 
      }
      
      void fill_halos_sgs_div(arr_t &a)
      {
        fill_halos_sclr(a);
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const int offset = 0)
      {
	using namespace idxperm;
        // the same logic as fill_halos_vctr_alng but have to consider offset ... TODO: find a way to reuse !
        if(!this->is_cyclic) 
        {
          if(halo == 1)
            this->recv(av[0 + offset],  idx_t(idx_ctor_arg_t(this->rght_halo_vctr))); 
          else
            this->xchng(av[0 + offset], idx_t(idx_ctor_arg_t(((this->rght_intr_vctr + off)^h)^(-1))), idx_t(idx_ctor_arg_t(this->rght_halo_vctr))); 
        }
        else
          this->xchng(av[0 + offset], idx_t(idx_ctor_arg_t(this->rght_intr_vctr + off)), idx_t(idx_ctor_arg_t(this->rght_halo_vctr))); 
      }
      
      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const real_t)
      {
        fill_halos_vctr_alng(av);
      }

      void fill_halos_vctr_nrml(arr_t &a)                 
      {                                                                         
        fill_halos_sclr(a);                                                  
      }  
      
      void fill_halos_vctr_alng_cyclic(arrvec_t<arr_t> &av, const bool ad = false)
      {
        fill_halos_vctr_alng(av, ad);
      }

      void fill_halos_vctr_nrml_cyclic(arr_t &a)
      {
        fill_halos_vctr_nrml(a);
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
