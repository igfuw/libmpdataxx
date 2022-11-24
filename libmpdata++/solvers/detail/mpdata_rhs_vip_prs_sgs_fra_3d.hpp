/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_fra_common.hpp>
#include <libmpdata++/formulae/fractal_reconstruction.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template<typename ct_params_t, int minhalo, class enableif = void>
      class mpdata_rhs_vip_prs_sgs_fra_dim
      {
//        static_assert(false);
      };

      template<typename ct_params_t, int minhalo>
      class mpdata_rhs_vip_prs_sgs_fra_dim<
        ct_params_t,
        minhalo,
        typename std::enable_if_t<ct_params_t::n_dims == 3>
      > : public detail::mpdata_rhs_vip_prs_sgs_fra_common<ct_params_t, minhalo>
      {
        public:

        using parent_t = detail::mpdata_rhs_vip_prs_sgs_fra_common<ct_params_t, minhalo>;
//        using parent_t::parent_t; // inheriting constructors
        using real_t = typename ct_params_t::real_t;
        using gen_t = typename std::conditional_t<std::is_same_v<real_t, float>, std::mt19937, std::mt19937_64>; // NOTE: we assume that if real_t is not float it has to be double

        struct ctor_args_t : parent_t::ctor_args_t
        {
          typename parent_t::bcp_ref_t
            &bcxl_ref, &bcxr_ref,
            &bcyl_ref, &bcyr_ref,
            &bczl_ref, &bczr_ref;
        };


        private:
  
        // helper variables for grid refinement
        // TODO: replace with idx_t: mid_ijk_r2r ijk_r2r_h
        rng_t mid_ijk_r2r_0, mid_ijk_r2r_1, mid_ijk_r2r_2;     // positions between already known values (to be filled during given iteration)
        rng_t ijk_r2r_0_h_with_halo, ijk_r2r_1_h, ijk_r2r_2_h; // all positions at resolution of given iteration
        int stride, hstride;                                   // stride and half stride

        // fill distmem halos of refinee, only at points overlapping with resolved points
        // and not outside of the modeled domain (e.g. in periodic case)
        // TODO: move to bcond or sth? would be filled only by remote bcond
        void fill_refinee_r2r_distmem_halos(const int e, const int halo_size)
        {
          const int e_ref = this->ix_r2r.at(e);
          // TODO: we only need to xchng along distmem direction (x)
          this->xchng(e);

          // TODO: assert halo_size <= solver halo size

          // TODO: DRY!
          switch(halo_size)
          {
            case 2:
              if(this->mem->distmem.rank() < this->mem->distmem.size() - 1)
                this->mem->psi_ref[e_ref](
                  this->mem->grid_size_ref[0].last() + this->mem->n_ref / 2 + this->mem->n_ref ,
                  this->ijk_r2r[1],
                  this->ijk_r2r[2]
                ) = 
                  this->mem->psi[e][0](
                    this->ijk[0].last()+2,
                    this->ijk[1],
                    this->ijk[2]
                  );
              if(this->mem->distmem.rank() > 0)
                this->mem->psi_ref[e_ref](
                  this->mem->grid_size_ref[0].first() - this->mem->n_ref / 2 - this->mem->n_ref ,
                  this->ijk_r2r[1],
                  this->ijk_r2r[2]
                ) = 
                  this->mem->psi[e][0](
                    this->ijk[0].first()-2,
                    this->ijk[1],
                    this->ijk[2]
                  );
              // no break intentionally
            case 1:
              if(this->mem->distmem.rank() < this->mem->distmem.size() - 1)
                this->mem->psi_ref[e_ref](
                  this->mem->grid_size_ref[0].last() + this->mem->n_ref / 2,
                  this->ijk_r2r[1],
                  this->ijk_r2r[2]
                ) = 
                  this->mem->psi[e][0](
                    this->ijk[0].last()+1,
                    this->ijk[1],
                    this->ijk[2]
                  );
              if(this->mem->distmem.rank() > 0)
                this->mem->psi_ref[e_ref](
                  this->mem->grid_size_ref[0].first() - this->mem->n_ref / 2,
                  this->ijk_r2r[1],
                  this->ijk_r2r[2]
                ) = 
                  this->mem->psi[e][0](
                    this->ijk[0].first()-1,
                    this->ijk[1],
                    this->ijk[2]
                  );
              break;
            default:
              assert(false);
              break;
          }
        }

        void refinement_ranges(const int iter, const int halo_size)
        {
          // messy, because in domain decomposition (sharedmem and distmem) some refined scalars are on the edge of the subdomain...
          if(iter==0)
          {
            mid_ijk_r2r_0 = this->rng_midpoints(this->ijk_r2r[0], this->mem->distmem.rank(), this->mem->distmem.size());
            mid_ijk_r2r_1 = this->rng_midpoints_in_out(this->ijk_r2r[1], this->rank, this->mem->size); // different because we dont want overlapping ranges in y 
            mid_ijk_r2r_2 = this->rng_midpoints(this->ijk_r2r[2]);

            ijk_r2r_1_h = this->ijk_r2r[1];
            ijk_r2r_2_h = this->ijk_r2r[2];
          }
          else
          {
            if(iter==1)
              mid_ijk_r2r_0 = this->rng_midpoints_in(mid_ijk_r2r_0, this->mem->distmem.rank(), this->mem->distmem.size());
            else
              mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0);

            mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1);
            mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);

            ijk_r2r_1_h = this->rng_half_stride(ijk_r2r_1_h, this->rank, this->mem->size);
            ijk_r2r_2_h = this->rng_half_stride(ijk_r2r_2_h);
          }

          stride = ijk_r2r_1_h.stride();
          assert(ijk_r2r_1_h.stride() == ijk_r2r_2_h.stride());
          assert(stride % 2 == 0);
          hstride = stride / 2;

          ijk_r2r_0_h_with_halo = rng_t(
            this->mem->distmem.rank() == 0                             ? this->ijk_r2r[0].first() : this->ijk_r2r[0].first() - halo_size * this->mem->n_ref,
            this->mem->distmem.rank() == this->mem->distmem.size() - 1 ? this->ijk_r2r[0].last()  : this->ijk_r2r[0].last()  + halo_size * this->mem->n_ref,
            hstride
          );
        }

        protected:

        // TODO: stretching parameters at the overlaps of the reconstructed and resolved arrays (ijk_r2r) are not used, do not generate them
        //       also not all parameters in the halo are needed (but some are!)
        void generate_stretching_parameters(const typename gen_t::result_type rng_seed)
        {
          gen_t gen(rng_seed); 
          std::uniform_real_distribution<> dis(-1, 1); // [-1,1), but whatever
          auto rand = std::bind(dis, gen);
          std::generate(this->d_j(this->ijk_ref_with_halo).begin(), this->d_j(this->ijk_ref_with_halo).end(), rand);
          this->d_j(this->ijk_ref_with_halo) = formulae::fractal::d_of_CDF_fctr<real_t>{}(this->d_j(this->ijk_ref_with_halo));
          this->mem->barrier();
        }

        // calculate refined points using (tri?)linear interpolation
        void interpolate_refinee(const int e = 0)
        {
          assert(opts::isset(ct_params_t::fractal_recon, opts::bit(e)));
          const int halo_size = 1;
          const int e_ref = this->ix_r2r.at(e);

          //// TEMPORARY
          //this->mem->psi_ref[e_ref] = -1000;
          //this->mem->barrier();


          fill_refinee_r2r_distmem_halos(e, halo_size);

          // fill refined array at position where it overlaps with the resolved array
          this->mem->refinee(e_ref)(this->ijk_r2r) = this->mem->advectee(e)(this->ijk);

//          generate_stretching_parameters();

          for(int i=0; i<this->n_fra_iter; ++i)
          {
            refinement_ranges(i, halo_size);

            formulae::fractal::intrp<0, real_t>(this->mem->psi_ref[e_ref], mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            this->mem->barrier();
            formulae::fractal::intrp<1, real_t>(this->mem->psi_ref[e_ref], mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h_with_halo, hstride);
            formulae::fractal::intrp<2, real_t>(this->mem->psi_ref[e_ref], mid_ijk_r2r_2, ijk_r2r_0_h_with_halo, this->rng_merge(ijk_r2r_1_h, mid_ijk_r2r_1), hstride);
          }
          this->mem->barrier();
        }

        void reconstruct_refinee(const int e = 0)
        {
          assert(opts::isset(ct_params_t::fractal_recon, opts::bit(e)));
          const int halo_size = 2;
          const int e_ref = this->ix_r2r.at(e);

          //// TEMPORARY
          //this->mem->psi_ref[e_ref] = -1000;
          //this->mem->barrier();

          //std::cerr << "mem->grid_size_ref[0]: " << this->mem->grid_size_ref[0] << std::endl;
          //std::cerr << "mem->grid_size_ref[1]: " << this->mem->grid_size_ref[1] << std::endl;
          //std::cerr << "mem->grid_size_ref[2]: " << this->mem->grid_size_ref[2] << std::endl;

          fill_refinee_r2r_distmem_halos(e, halo_size);

          // fill refined array at position where it overlaps with the resolved array
          this->mem->refinee(e_ref)(this->ijk_r2r) = this->mem->advectee(e)(this->ijk);

//          generate_stretching_parameters();
//          this->mem->barrier();

          for(int i=0; i<this->n_fra_iter; ++i)
          {
            refinement_ranges(i, halo_size);

            formulae::fractal::rcnstrct<0, real_t>(this->mem->psi_ref[e_ref], this->rng_dbl_stride(mid_ijk_r2r_0), ijk_r2r_1_h,           ijk_r2r_2_h,                                 hstride, this->c_j, this->d_j, this->f_j);
            this->mem->barrier();
            formulae::fractal::rcnstrct<1, real_t>(this->mem->psi_ref[e_ref], this->rng_dbl_stride(mid_ijk_r2r_1)        , ijk_r2r_2_h,           ijk_r2r_0_h_with_halo,                       hstride, this->c_j, this->d_j, this->f_j); // NOTE: rng_dbl_stride(mid_ijk_r2r_1) gives overlapping ranges between thread subdomains... however it seems to work and naive fixes didnt work
            formulae::fractal::rcnstrct<2, real_t>(this->mem->psi_ref[e_ref], this->rng_dbl_stride(mid_ijk_r2r_2)        , ijk_r2r_0_h_with_halo, this->rng_merge(ijk_r2r_1_h, mid_ijk_r2r_1), hstride, this->c_j, this->d_j, this->f_j);
          }
          this->mem->barrier();
        }

        void xchng_ref(
                       typename parent_t::arr_t &arr,
                       const idx_t<3> &range_ijk
        )
        {
          this->mem->barrier();
          for (auto &bc : this->bcs_ref[1]) bc->fill_halos_sclr(arr, range_ijk[2], range_ijk[0]);
          for (auto &bc : this->bcs_ref[2]) bc->fill_halos_sclr(arr, range_ijk[0], range_ijk[1]);
          for (auto &bc : this->bcs_ref[0]) bc->fill_halos_sclr(arr, range_ijk[1], range_ijk[2]);
          this->mem->barrier();
        }


        void xchng_vctr_ref(
                       arrvec_t<typename parent_t::arr_t> &arrvec,
                       const idx_t<3> &range_ijk
        )
        {
          this->mem->barrier();
          for (auto &bc : this->bcs_ref[1]) bc->fill_halos_sclr(arrvec[1], range_ijk[2], range_ijk[0]);
          for (auto &bc : this->bcs_ref[2]) bc->fill_halos_sclr(arrvec[2], range_ijk[0], range_ijk[1]);
          for (auto &bc : this->bcs_ref[0]) bc->fill_halos_sclr(arrvec[0], range_ijk[1], range_ijk[2]);
          this->mem->barrier();
        }

//        void xchng_ref(int e)
//        {
//          const int e_ref = this->ix_r2r.at(e);
//          this->xchng_sclr_ref(this->mem->psi_ref[e_ref], this->ijk_ref);
//        }

        // calculate courant numbers at refined cell edges using a
        // linear interpolation of a refined uvw (src) from cell centers to cell edges (dste
        void interpolate_refined_courants(arrvec_t<typename parent_t::arr_t> &dst,
                                          arrvec_t<typename parent_t::arr_t> &src)
        {
          this->xchng_vctr_ref(src, this->ijk_ref);
          this->template intrp<0>(dst[0], src[0], this->ijk_ref[0]^h, this->ijk_ref[1], this->ijk_ref[2], this->dijk_ref[0], true);
          this->template intrp<1>(dst[1], src[1], this->ijk_ref[1]^h, this->ijk_ref[2], this->ijk_ref[0], this->dijk_ref[1], true);
          this->template intrp<2>(dst[2], src[2], this->ijk_ref[2]^h, this->ijk_ref[0], this->ijk_ref[1], this->dijk_ref[2], true);
          this->mem->barrier();
        }

        public:

        mpdata_rhs_vip_prs_sgs_fra_dim(
          ctor_args_t args,
          const typename parent_t::rt_params_t &p
        ) :
          parent_t(args, p)
        {
          this->dijk_ref = {this->di / this->n_ref, this->dj / this->n_ref, this->dk / this->n_ref};
          this->ijkm_ref = this->ijk_ref;
          this->ijkm_ref.lbound(0) = this->ijkm_ref.lbound(0) - 1;
          this->ijkm_ref.lbound(1) = this->ijkm_ref.lbound(1) - 1;
          this->ijkm_ref.lbound(2) = this->ijkm_ref.lbound(2) - 1;

          this->set_bcs(this->bcs_ref, 0, args.bcxl_ref, args.bcxr_ref);
          this->set_bcs(this->bcs_ref, 1, args.bcyl_ref, args.bcyr_ref);
          this->set_bcs(this->bcs_ref, 2, args.bczl_ref, args.bczr_ref);

        /*
          this->bcs_ref[0][0].reset(
            new bcond::bcond<real_t, parent_t::halo_ref, bcond::open, bcond::left, ct_params_t::n_dims, 0>(
              this->ijk_ref[0],
              this->mem->distmem.grid_size_ref
            )
          );
          */

        /*
          this->set_bcs(this->bcs_ref, 0, args.bcxl_ref, args.bcxr_ref);
          this->set_bcs(this->bcs_ref, 1, args.bcyl_ref, args.bcyr_ref);
          this->set_bcs(this->bcs_ref, 2, args.bczl_ref, args.bczr_ref);
          */
        }

        // helper method to allocate n_arr refined scalar temporary arrays
        static void alloc_tmp_sclr_ref(
          typename parent_t::mem_t *mem,
          const char * __file__, const int n_arr,
          std::string name = "",
          bool srfc = false
        )
        {
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());

          if (!name.empty()) mem->avail_tmp[name] = std::make_pair(__file__, mem->tmp[__file__].size() - 1);

          for (int n = 0; n < n_arr; ++n)
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
              parent_t::rng_ref_distmem_halo(mem->grid_size_ref[0]^parent_t::halo_ref, mem->n_ref, mem->distmem.rank(), mem->distmem.size()), 
              mem->grid_size_ref[1]^parent_t::halo_ref,
              srfc ? rng_t(0, 0) : mem->grid_size_ref[2]^parent_t::halo_ref,
              arr3D_storage
            )));
        }

        static void alloc_tmp_vctr_ref(
          typename parent_t::mem_t *mem,
          const char * __file__
        )
        {
          const int n_arr = 3;
          const std::vector<std::vector<bool>> stgr{{true, false, false}, {false, true, false}, {false, false, true}};

          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          for (int n = 0; n < n_arr; ++n)
          {
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
              stgr[n][0] ? mem->grid_size_ref[0]^h : mem->grid_size_ref[0],
              stgr[n][1] ? mem->grid_size_ref[1]^h : mem->grid_size_ref[1],
              stgr[n][2] ? mem->grid_size_ref[2]^h : mem->grid_size_ref[2],
              arr3D_storage
            )));
          }
        }

      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
