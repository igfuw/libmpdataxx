// 2D MPI ``remote'' boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/remote_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>    
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == remote && 
        dir == left   && 
        n_dims == 2
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {
      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : -1;

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
        using namespace idxperm;
        this->xchng(a, pi<d>(this->left_intr_sclr + off, j), pi<d>(this->left_halo_sclr, j));
      }

      void fill_halos_pres(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }
      
      void save_edge_vel(const arr_t &, const rng_t &) {}

      void set_edge_pres(arr_t &, const rng_t &, int) {}

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
        using namespace idxperm;
        if(!this->is_cyclic) 
        {
          if(halo == 1)
            // send vectors to the left of the domain
            this->send(av[0], pi<d>(this->left_intr_vctr + off, j));
          else
            // receive the halo without the rightmost column, which was caluclated by this process
            this->xchng(av[0], pi<d>(this->left_intr_vctr + off, j), pi<d>((this->left_halo_vctr^h)^(-1), j));
        }
        else
          this->xchng(av[0], pi<d>(this->left_intr_vctr + off, j), pi<d>(this->left_halo_vctr, j));
      }

      // TODO: sgs fill_halos

      // TODO: move to common? (same in cyclic!)
      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)                 
      {                                                                         
        fill_halos_sclr(a, j);                                                  
      }  
    };

    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == remote &&
        dir == rght   &&
        n_dims == 2
      >::type
    > : public detail::remote_common<real_t, halo, dir, n_dims>
    {
      using parent_t = detail::remote_common<real_t, halo, dir, n_dims>;
      using arr_t = blitz::Array<real_t, 2>;
      using parent_t::parent_t; // inheriting ctor

      const int off = this->is_cyclic ? 0 : 1;

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const bool deriv = false)
      {
        using namespace idxperm;
        this->xchng(a, pi<d>(this->rght_intr_sclr + off, j), pi<d>(this->rght_halo_sclr, j));
      }

      void fill_halos_pres(arr_t &a, const rng_t &j)
      {
        fill_halos_sclr(a, j);
      }
      
      void save_edge_vel(const arr_t &, const rng_t &) {}

      void set_edge_pres(arr_t &, const rng_t &, int) {}

      void fill_halos_vctr_alng(arrvec_t<arr_t> &av, const rng_t &j, const bool ad = false)
      {
        using namespace idxperm;
        if(!this->is_cyclic)
        {
          if(halo == 1)
            //receive the halo
            this->recv(av[0], pi<d>(this->rght_halo_vctr, j));
          else
            // don't send the first column to the right of the domain, it will be calculated and sent here by the process to the right
            this->xchng(av[0], pi<d>(((this->rght_intr_vctr + off)^h)^(-1), j), pi<d>(this->rght_halo_vctr, j));
        }
        else
          this->xchng(av[0], pi<d>(this->rght_intr_vctr + off, j), pi<d>(this->rght_halo_vctr, j));
      }

      // TODO: move to common? (same in cyclic!)
      void fill_halos_vctr_nrml(arr_t &a, const rng_t &j)                 
      {                                                                         
        fill_halos_sclr(a, j);                                                  
      }  
    };
  } // namespace bcond
} // namespace libmpdataxx
