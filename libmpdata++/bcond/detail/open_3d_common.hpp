// common code for ``open'' MPI boundary conditions for libmpdata++
//
// licensing: GPU GPL v3
// copyright: University of Warsaw
//
//#pragma once
//
//#include <libmpdata++/bcond/detail/bcond_common.hpp>
//
//namespace libmpdataxx
//{
//  namespace bcond
//  {
//    namespace detail
//    {
//      template <typename real_t, int halo>
//      class open_3d_common : public bcond_common<real_t, halo, 3>
//      {
//        using parent_t = detail::bcond_common<real_t, halo, 3>;
//
//        protected: 
//
//        const int thread_rank, thread_size;
//
//        public:
//
//        // ctor
//        open_3d_common(
//          const rng_t &i,
//          const std::array<int, 3> &distmem_grid_size,
//          const int thread_rank,
//          const int thread_size
//        ) :
//          parent_t(i, distmem_grid_size),
//          thread_rank(thread_rank),
//          thread_size(thread_size)
//        {}
//
//        open_3d_common(
//          const rng_t &i,
//          const std::array<int, 3> &distmem_grid_size
//        ) :
//          parent_t(i, distmem_grid_size),
//          thread_rank(0),
//          thread_size(0)
//        {
//          throw std::runtime_error("Wrong open_3d_common ctor called.");
//        }
//      };
//    };
//  } // namespace bcond
//} // namespace libmpdataxx
