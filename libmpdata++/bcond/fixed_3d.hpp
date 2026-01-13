// 3D fixed (Dirichlet) boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/bcond/detail/bcond_common.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == fixed &&
        dir == left &&
        n_dims == 3
      >::type
    > : public bcond<real_t, halo, open, dir, n_dims, d> // TODO: in some cases (e.g. cloud chamber) we might want to inherit from rigid?
    {
      using parent_t = bcond<real_t, halo, open, dir, n_dims, d>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      // holds fixed wall values (equal to initial edge values)
      std::unordered_map<const arr_t*, arr_t> edge_values;

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        #if !defined(NDEBUG)
          assert(edge_values.find(&a) != edge_values.end() && "fixed bcond: edge values not saved before filling halos");
        #endif

        using namespace idxperm;
        for (int i = this->left_halo_sclr.first(); i <= this->left_halo_sclr.last(); ++i)
        {
          a(pi<d>(i, j, k)) = edge_values[&a](pi<d>(0, j, k));
        }
      }

      void save_edge_val(const arr_t &a, const arr_t &val, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // assert(a.shape() == val.shape());
        auto s = a.shape();
        s[d] = 1;
        edge_values.emplace(&a, arr_t(s));
        if constexpr (d == 0)      edge_values[&a].reindexSelf({0, a.lbound(1), a.lbound(2)});
        else if constexpr (d == 1) edge_values[&a].reindexSelf({a.lbound(0), 0, a.lbound(2)});
        else if constexpr (d == 2) edge_values[&a].reindexSelf({a.lbound(0), a.lbound(1), 0});
        edge_values[&a](pi<d>(0, j, k)) = val(pi<d>(this->left_edge_sclr, j, k));
      }
    };


    template <typename real_t, int halo, bcond_e knd, drctn_e dir, int n_dims, int d>
    class bcond<       real_t,     halo,         knd,         dir,     n_dims,     d,
      typename std::enable_if<
        knd == fixed &&
        dir == rght &&
        n_dims == 3
      >::type
    > : public bcond<real_t, halo, open, dir, n_dims, d> // TODO: in some cases (e.g. cloud chamber) we might want to inherit from rigid?
    {
      using parent_t = bcond<real_t, halo, open, dir, n_dims, d>;
      using arr_t = blitz::Array<real_t, 3>;
      using parent_t::parent_t; // inheriting ctor

      // holds fixed wall values (equal to initial edge values)
      std::unordered_map<const arr_t*, arr_t> edge_values;

      public:

      void fill_halos_sclr(arr_t &a, const rng_t &j, const rng_t &k, const bool deriv = false)
      {
        #if !defined(NDEBUG)
          assert(edge_values.find(&a) != edge_values.end() && "fixed bcond: edge values not saved before filling halos");
        #endif

        using namespace idxperm;
        for (int i = this->rght_halo_sclr.first(); i <= this->rght_halo_sclr.last(); ++i)
        {
          a(pi<d>(i, j, k)) = edge_values[&a](pi<d>(0, j, k));
        }
      }

      void save_edge_val(const arr_t &a, const arr_t &val, const rng_t &j, const rng_t &k)
      {
        using namespace idxperm;
        // assert(a.shape() == val.shape());
        auto s = a.shape();
        s[d] = 1;
        edge_values.emplace(&a, arr_t(s));
        if constexpr (d == 0)      edge_values[&a].reindexSelf({0, a.lbound(1), a.lbound(2)});
        else if constexpr (d == 1) edge_values[&a].reindexSelf({a.lbound(0), 0, a.lbound(2)});
        else if constexpr (d == 2) edge_values[&a].reindexSelf({a.lbound(0), a.lbound(1), 0});
        edge_values[&a](pi<d>(0, j, k)) = val(pi<d>(this->rght_edge_sclr, j, k));
      }
    };
  } // namespace bcond
} // namespace libmpdataxx
