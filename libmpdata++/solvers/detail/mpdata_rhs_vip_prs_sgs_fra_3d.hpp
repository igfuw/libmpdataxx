/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_fra_common.hpp>

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
        using parent_t = detail::mpdata_rhs_vip_prs_sgs_fra_common<ct_params_t, minhalo>;
        using parent_t::parent_t; // inheriting constructors

        protected:

        // interpolation similar to mpdata_rhs_vip...
        template<int d, class arr_t>
        void intrp(
          arr_t arr,
          const rng_t &i,
          const rng_t &j,
          const rng_t &k,
          const int &dist
        )
        {
          using idxperm::pis;
//          using namespace arakawa_c;
          using real_t = typename ct_params_t::real_t;

//          if(d==0)
//            std::cerr << "range<" << d << ">: " << i << " " << j << " " << k << std::endl;
//          if(d==1)
//            std::cerr << "range<" << d << ">: " << k << " " << i << " " << j << std::endl;
//          if(d==2)
//            std::cerr << "range<" << d << ">: " << j << " " << k << " " << i << std::endl;

//          std::cerr << "range - dist: " << i - dist << " " << j << " " << k << std::endl;
//          std::cerr << "range + dist: " << i + dist << " " << j << " " << k << std::endl;
//
//          std::cerr << "arr range: " << arr(i, j, k) << std::endl;
//          std::cerr << "arr range - dist: " << arr(i - dist, j, k) << std::endl;
//          std::cerr << "arr range + dist: " << arr(i + dist, j, k) << std::endl;
//
//          std::cerr << "arr range pi: " << arr(pis<d>(i, j, k)) << std::endl;
//          std::cerr << "arr range - dist pi: " << arr(pis<d>(i - dist, j, k)) << std::endl;
//          std::cerr << "arr range + dist pi: " << arr(pis<d>(i + dist, j, k)) << std::endl;

          arr(pis<d>(i, j, k)) = real_t(.5) * (
            arr(pis<d>(i - dist, j, k)) +
            arr(pis<d>(i + dist, j, k))
          );
        }

        void interpolate_refinee(const int e = 0)
        {
          using namespace arakawa_c; // for rng_t operator^

          rng_t mid_ijk_r2r_0, mid_ijk_r2r_1, mid_ijk_r2r_2; // positions between already known values (to be filled during given iteration)
          rng_t ijk_r2r_0_h, ijk_r2r_1_h, ijk_r2r_2_h;       // all positions at resolution of given iteration
          int stride, hstride;

          // fill refined array at position where it overlaps with the resolved array
          this->mem->refinee(e)(this->ijk_r2r) = this->mem->advectee(e)(this->ijk);

//          interpolate_refinee_on_edges(); // with MPI, some refined points are at the edges of the domain; calculate them using halos of non-refined arrays
// TODO: fill halo sclr ?

          for(int i=0; i<this->n_fra_iter; ++i)
          {
            // messy, because in domain decomposition (sharedmem and distmem) some refined scalars are on the edge of the subdomain...
            if(i==0)
            {
              mid_ijk_r2r_0 = this->rng_midpoints(this->ijk_r2r[0], this->mem->distmem.rank(), this->mem->distmem.size());
              mid_ijk_r2r_1 = this->rng_midpoints(this->ijk_r2r[1], this->rank, this->mem->size);
              mid_ijk_r2r_2 = this->rng_midpoints(this->ijk_r2r[2]);

              ijk_r2r_0_h = this->ijk_r2r[0];
              ijk_r2r_1_h = this->ijk_r2r[1];
              ijk_r2r_2_h = this->ijk_r2r[2];
            }
            else
            {
              if(i==1)
              {
                mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0, this->mem->distmem.rank(), this->mem->distmem.size());
                mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1, this->rank, this->mem->size);
                mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);
              }
              else
              {
                mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0);
                mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1);
                mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);
              }

              ijk_r2r_0_h = this->rng_half_stride(ijk_r2r_0_h);
              ijk_r2r_1_h = this->rng_half_stride(ijk_r2r_1_h);
              ijk_r2r_2_h = this->rng_half_stride(ijk_r2r_2_h);
            }

            stride = ijk_r2r_0_h.stride();
            assert(ijk_r2r_0_h.stride() == ijk_r2r_1_h.stride() == ijk_r2r_2_h.stride());
            assert(stride % 2 == 0);
            hstride = stride / 2;

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, ijk_r2r_1_h, hstride);

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, mid_ijk_r2r_1, ijk_r2r_2_h, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, ijk_r2r_1_h, hstride);

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, mid_ijk_r2r_2, hstride);
//            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, mid_ijk_r2r_0, hstride);
//            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, mid_ijk_r2r_1, hstride);

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, mid_ijk_r2r_1, mid_ijk_r2r_2, hstride);
//            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, mid_ijk_r2r_0, hstride);
//            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);
          }
        }

        public:

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
              mem->grid_size_ref[0], // NOTE: no halo
              mem->grid_size_ref[1],
              srfc ? rng_t(0, 0) : mem->grid_size_ref[2],
              arr3D_storage
            )));
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
