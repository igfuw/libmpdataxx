// 3D ground/sky boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    // ground
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == gndsky &&
        dir == left &&
        n_dims == 3
      >::type
    > : public bcond<real_t, halo, rigid, dir, n_dims, d>
    { 
      using parent_t = bcond<real_t, halo, rigid, dir, n_dims, d>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      public:

      void fill_halos_sgs_div(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        a(pi<d>(this->left_edge_sclr - h, j, k)) = 2 * a(pi<d>(this->left_edge_sclr + h, j, k))
                                                     - a(pi<d>(this->left_edge_sclr + 1 + h, j, k));
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &b, const rng_t &j, const rng_t &k, const int offset = 0)
      {
        using namespace idxperm;
        // fill halos for a staggered field based on prescribed values on the ground (b array)
        // that is 0.5 * (a(gnd-h) + a(gnd+h)) = b(i)
        const auto &a = av[offset + d];
        a(pi<d>(this->left_edge_sclr - h, j, k)) = 2 * b(pi<d>(this->left_edge_sclr, j, k)) - a(pi<d>(this->left_edge_sclr + h, j, k));
      }
      
      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &u, const arr_t &div, const rng_t &j, const rng_t &k, const real_t di)
      {
        using namespace idxperm;
        const auto &a = av[d];
        a(pi<d>(this->left_edge_sclr - h, j, k)) = 2 * ( ( 3 * u(pi<d>(this->left_edge_sclr + 1, j, k))
                                                         - 2 * u(pi<d>(this->left_edge_sclr, j, k)) 
                                                         -     u(pi<d>(this->left_edge_sclr + 2, j, k)) 
                                                         ) / di
                                                       - div(pi<d>(this->left_edge_sclr - h, j, k))
                                                       );
      }
    };

    // sky
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,  
      typename std::enable_if<
        knd == gndsky &&
        dir == rght &&
        n_dims == 3
      >::type
    > : public bcond<real_t, halo, rigid, dir, n_dims, d>
    {
      using parent_t = bcond<real_t, halo, rigid, dir, n_dims, d>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor
      
      public:

      void fill_halos_sgs_div(arr_t &a, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        a(pi<d>(this->rght_edge_sclr + h, j, k)) = 2 * a(pi<d>(this->rght_edge_sclr - h, j, k))
                                                   -   a(pi<d>(this->rght_edge_sclr - 1 - h, j, k));
      }

      void fill_halos_sgs_vctr(arrvec_t<arr_t> &av, const arr_t &, const rng_t &j, const rng_t &k, const int offset = 0)
      {
        // fill halos for a staggered field so that it has zero value on tke sky edge
        // that is 0.5 * (a(sky-h) + a(sky+h)) = 0
        using namespace idxperm;
        const auto &a = av[offset + d];
        a(pi<d>(this->rght_edge_sclr + h, j, k)) = - a(pi<d>(this->rght_edge_sclr - h, j, k));
      }
      
      void fill_halos_sgs_tnsr(arrvec_t<arr_t> &av, const arr_t &, const arr_t &, const rng_t &j, const rng_t &k, const real_t di)
      {
        using namespace idxperm;
        const auto &a = av[d];
        a(pi<d>(this->rght_edge_sclr + h, j, k)) = 2 * a(pi<d>(this->rght_edge_sclr - h, j, k))
                                                   -   a(pi<d>(this->rght_edge_sclr - 1 - h, j, k));
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
