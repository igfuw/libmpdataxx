/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// for solvers with both a refined and regular grid

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
        bcond::bcond_e bczl, bcond::bcond_e bczr
      >
      class concurr_common_ref : public concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
      {
        private:

        using parent_t = concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
//        using parent_t::parent_t;

        boost::mpi::communicator mpic, mpic_ref;

        public:
        // ctor
        concurr_common_ref(
          const typename solver_t_::rt_params_t &p,
          typename solver_t_::mem_t *mem_p,
          const int &size
        ) :
          parent_t(p, mem_p, size),
          mpic(MPI_COMM_WORLD, boost::mpi::comm_duplicate),
          mpic_ref(MPI_COMM_WORLD, boost::mpi::comm_duplicate)
        {}

        protected:
        // 3D version, note sharedmem in y direction!
        void init(
          const typename solver_t_::rt_params_t &p,
          const int &n1, const int &n0 = 1, const int &n2 = 1
        ) {

          typename solver_t_::bcp_ref_t bxl_ref, bxr_ref, byl_ref, byr_ref, bzl_ref, bzr_ref, shrdl_ref, shrdr_ref;

          // TODO: renew pointers only if invalid ?
          for (int i0 = 0; i0 < n0; ++i0)
          {
            for (int i1 = 0; i1 < n1; ++i1)
            {
              for (int i2 = 0; i2 < n2; ++i2)
              {
                // i1 is the local thread rank, n1 is the number of threads. These are needed by remote bcond, because only rank=0 does mpi communication
                this->template init_bcs_3d<solver_t_::halo>(
                  this->bxl, this->bxr, this->byl, this->byr, this->bzl, this->bzr, this->shrdl, this->shrdr,
                  this->mem->grid_size,
                  this->mem->distmem.grid_size,
                  i1, n1, mpic);

                this->template init_bcs_3d<solver_t_::halo_ref>(
                  bxl_ref, bxr_ref, byl_ref, byr_ref, bzl_ref, bzr_ref, shrdl_ref, shrdr_ref,
                  this->mem->grid_size_ref,
                  this->mem->distmem.grid_size_ref,
                  i1, n1, mpic_ref);

                /*
                this->init_bcs_3d(i1, n1);

                this->template bc_set<bcxl, bcond::left, 0, solver_t_::halo_ref>(bxl_ref, this->mem->grid_size_ref, this->mem->distmem.grid_size_ref, i1, n1);
                this->template bc_set<bcxr, bcond::rght, 0, solver_t_::halo_ref>(bxr_ref, this->mem->grid_size_ref, this->mem->distmem.grid_size_ref, i1, n1);

                this->template bc_set<bcyl, bcond::left, 1, solver_t_::halo_ref>(byl_ref, this->mem->grid_size_ref, this->mem->distmem.grid_size_ref);
                this->template bc_set<bcyr, bcond::rght, 1, solver_t_::halo_ref>(byr_ref, this->mem->grid_size_ref, this->mem->distmem.grid_size_ref);

                this->template bc_set<bczl, bcond::left, 2, solver_t_::halo_ref>(bzl_ref, this->mem->grid_size_ref, this->mem->distmem.grid_size_ref);
                this->template bc_set<bczr, bcond::rght, 2, solver_t_::halo_ref>(bzr_ref, this->mem->grid_size_ref, this->mem->distmem.grid_size_ref);

                shrdl_ref.reset(new bcond::shared<typename solver_t_::real_t, solver_t_::halo_ref, solver_t_::n_dims>()); // TODO: shrdy if n1 != 1
                shrdr_ref.reset(new bcond::shared<typename solver_t_::real_t, solver_t_::halo_ref, solver_t_::n_dims>()); // TODO: shrdy if n1 != 1
                */

                this->algos.push_back(
                  new solver_t_(
                    typename solver_t_::ctor_args_t{
                      {
                        i1,
                        this->mem.get(),
                        this->bxl, this->bxr,
                        i1 == 0      ? this->byl : this->shrdl,
                        i1 == n1 - 1 ? this->byr : this->shrdr,
                        this->bzl, this->bzr,
                        this->mem->slab(this->mem->grid_size[0], i0, n0),
                        this->mem->slab(this->mem->grid_size[1], i1, n1),
                        this->mem->slab(this->mem->grid_size[2], i2, n2)
                      },
                      bxl_ref, bxr_ref,
                      i1 == 0      ? byl_ref : shrdl_ref,
                      i1 == n1 - 1 ? byr_ref : shrdr_ref,
                      bzl_ref, bzr_ref,
                    },
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
