/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// for solver with regular grid only
// needed because regular and refined grid solvers have different ctors (ref needs bconds for ref grid)

#pragma once

#include <libmpdata++/concurr/detail/concurr_common_hlpr.hpp>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr,
	class enable = void
      >
      class concurr_common_reg
      {};

      // 1D version
      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
      >
      class concurr_common_reg<
        solver_t_,
        bcxl, bcxr,
        bcyl, bcyr,
        bczl, bczr,
        typename std::enable_if_t<solver_t_::n_dims == 1>
      > : public concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
      {
        private:

        using parent_t = concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
        using parent_t::parent_t;
	using real_t = typename solver_t_::real_t;

        protected:

        void init(
          const typename solver_t_::rt_params_t &p,
          const int &n0
        )
        {
          // NOTE: for remote bcond, thread_rank set to 0 on purpose in 1D to have propre left/right message tags 
          this->template bc_set<bcxl, bcond::left, 0, solver_t_::halo>(this->bxl, this->mem->grid_size, this->mem->distmem.grid_size);
          this->template bc_set<bcxr, bcond::rght, 0, solver_t_::halo>(this->bxr, this->mem->grid_size, this->mem->distmem.grid_size);

          for (int i0 = 0; i0 < n0; ++i0)
          {
            this->shrdl.reset(new bcond::shared<real_t, solver_t_::halo, solver_t_::n_dims>());
            this->shrdr.reset(new bcond::shared<real_t, solver_t_::halo, solver_t_::n_dims>());

            this->algos.push_back(
              new solver_t_(
                typename solver_t_::ctor_args_t({
                  i0,
                  this->mem.get(),
                  i0 == 0      ? this->bxl : this->shrdl,
                  i0 == n0 - 1 ? this->bxr : this->shrdr,
                  this->mem->slab(this->mem->grid_size[0], i0, n0)
                }),
                p
              )
            );
          }
        }
      };


      // 2D version
      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
      >
      class concurr_common_reg<
        solver_t_,
        bcxl, bcxr,
        bcyl, bcyr,
        bczl, bczr,
        typename std::enable_if_t<solver_t_::n_dims == 2>
      > : public concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
      {
        private:

        using parent_t = concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
        using parent_t::parent_t;
	using real_t = typename solver_t_::real_t;

        protected:
        // TODO: assert parallelisation in the right dimensions! (blitz::assertContiguous)
        void init(
          const typename solver_t_::rt_params_t &p,
          const int &n0, const int &n1 = 1
        ) {
          for (int i0 = 0; i0 < n0; ++i0)
          {
            for (int i1 = 0; i1 < n1; ++i1)
            {
              // NOTE: for remote bcond, thread_rank set to 0 on purpose in 2D to have propre left/right message tags 
              this->template bc_set<bcxl, bcond::left, 0, solver_t_::halo>(this->bxl, this->mem->grid_size, this->mem->distmem.grid_size);
              this->template bc_set<bcxr, bcond::rght, 0, solver_t_::halo>(this->bxr, this->mem->grid_size, this->mem->distmem.grid_size);

              this->template bc_set<bcyl, bcond::left, 1, solver_t_::halo>(this->byl, this->mem->grid_size, this->mem->distmem.grid_size);
              this->template bc_set<bcyr, bcond::rght, 1, solver_t_::halo>(this->byr, this->mem->grid_size, this->mem->distmem.grid_size);

              this->shrdl.reset(new bcond::shared<real_t, solver_t_::halo, solver_t_::n_dims>()); // TODO: shrdy if n1 != 1
              this->shrdr.reset(new bcond::shared<real_t, solver_t_::halo, solver_t_::n_dims>()); // TODO: shrdy if n1 != 1

              this->algos.push_back(
                new solver_t_(
                  typename solver_t_::ctor_args_t({
                    i0,
                    this->mem.get(),
                    i0 == 0      ? this->bxl : this->shrdl,
                    i0 == n0 - 1 ? this->bxr : this->shrdr,
                    this->byl, this->byr,
                    this->mem->slab(this->mem->grid_size[0], i0, n0),
                    this->mem->slab(this->mem->grid_size[1], i1, n1)
                  }),
                  p
                )
              );
            }
          }
        }
      };

      // 3D version, note sharedmem in y direction!
      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
      >
      class concurr_common_reg<
        solver_t_,
        bcxl, bcxr,
        bcyl, bcyr,
        bczl, bczr,
        typename std::enable_if_t<solver_t_::n_dims == 3>
      > : public concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
      {
        private:

        using parent_t = concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
        using parent_t::parent_t;
	using real_t = typename solver_t_::real_t;

        protected:
        void init(
          const typename solver_t_::rt_params_t &p,
          const int &n1, const int &n0 = 1, const int &n2 = 1
        ) {
          // TODO: renew pointers only if invalid ?
          for (int i0 = 0; i0 < n0; ++i0)
          {
            for (int i1 = 0; i1 < n1; ++i1)
            {
              for (int i2 = 0; i2 < n2; ++i2)
              {
                this->template init_bcs_3d<solver_t_::halo>(
                  this->bxl, this->bxr, this->byl, this->byr, this->bzl, this->bzr, this->shrdl, this->shrdr, 
                  this->mem->grid_size, 
                  this->mem->distmem.grid_size, 
                  i1, n1);

                this->algos.push_back(
                  new solver_t_(
                    typename solver_t_::ctor_args_t({
                      i1,
                      this->mem.get(),
                      this->bxl, this->bxr,
                      i1 == 0      ? this->byl : this->shrdl,
                      i1 == n1 - 1 ? this->byr : this->shrdr,
                      this->bzl, this->bzr,
                      this->mem->slab(this->mem->grid_size[0], i0, n0),
                      this->mem->slab(this->mem->grid_size[1], i1, n1),
                      this->mem->slab(this->mem->grid_size[2], i2, n2)
                    }),
                    p
                  )
                );
              }
            }
          }
        }
      };
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
