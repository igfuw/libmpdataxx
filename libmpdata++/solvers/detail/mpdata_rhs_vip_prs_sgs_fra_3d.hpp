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

//          std::cerr << "range: " << i << " " << j << " " << k << std::endl;
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
//        arrvec_t<typename parent_t::arr_t> &dst,
//                                 const arrvec_t<typename parent_t::arr_t> &src,
//                                 const int stride)
        {
          using namespace arakawa_c; // for rng_t operator^

//          interpolate_refinee_on_edges(); // with MPI, some refined points are at the edges of the domain; calculate them using halos of non-refined arrays

          // first round of interpolation
          auto stride = this->ijk_r2r[0].stride();
          assert(stride % 2 == 0);
          assert(this->ijk_r2r[0].stride() == this->ijk_r2r[1].stride() == this->ijk_r2r[2].stride());
          auto hstride = stride / 2;

          const auto mid_ijk_r2r_0 = this->rng_midpoints(this->ijk_r2r[0]);
          const auto mid_ijk_r2r_1 = this->rng_midpoints(this->ijk_r2r[1]);
          const auto mid_ijk_r2r_2 = this->rng_midpoints(this->ijk_r2r[2]);

          // interpolate between points of the large grid
          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, this->ijk_r2r[1], this->ijk_r2r[2], hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, this->ijk_r2r[2], this->ijk_r2r[0], hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, this->ijk_r2r[0], this->ijk_r2r[1], hstride);

          // interpolate between just interpolated points
          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, mid_ijk_r2r_1, this->ijk_r2r[2], hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, this->ijk_r2r[0], hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, this->ijk_r2r[1], hstride);

          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, this->ijk_r2r[1], mid_ijk_r2r_2, hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, this->ijk_r2r[2], mid_ijk_r2r_0, hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, this->ijk_r2r[0], mid_ijk_r2r_1, hstride);

          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, mid_ijk_r2r_1, mid_ijk_r2r_2, hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, mid_ijk_r2r_0, hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);

          // subsequent rounds of interpolation ...

          stride = mid_ijk_r2r_0_2.stride();
          assert(stride % 2 == 0);
          assert(mid_ijk_r2r_0_2.stride() == mid_ijk_r2r_1_2.stride() == mid_ijk_r2r_2_2.stride());
          hstride = stride / 2;

          const auto mid_ijk_r2r_0_2 = this->rng_midpoints_out(mid_ijk_r2r_0);
          const auto mid_ijk_r2r_1_2 = this->rng_midpoints_out(mid_ijk_r2r_1);
          const auto mid_ijk_r2r_2_2 = this->rng_midpoints_out(mid_ijk_r2r_2);

          const auto ijk_r2r_0_h = this->rng_half_stride(this->ijk_r2r[0]);
          const auto ijk_r2r_1_h = this->rng_half_stride(this->ijk_r2r[1]);
          const auto ijk_r2r_2_h = this->rng_half_stride(this->ijk_r2r[2]);

          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0_2, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1_2, ijk_r2r_2_h, ijk_r2r_0_h, hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2_2, ijk_r2r_0_h, ijk_r2r_1_h, hstride);

          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0_2, mid_ijk_r2r_1_2, ijk_r2r_2_h, hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1_2, mid_ijk_r2r_2_2, ijk_r2r_0_h, hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2_2, mid_ijk_r2r_0_2, ijk_r2r_1_h, hstride);

          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0_2, ijk_r2r_1_h, mid_ijk_r2r_2_2, hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1_2, ijk_r2r_2_h, mid_ijk_r2r_0_2, hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2_2, ijk_r2r_0_h, mid_ijk_r2r_1_2, hstride);

          intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0_2, mid_ijk_r2r_1_2, mid_ijk_r2r_2_2, hstride);
          intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1_2, mid_ijk_r2r_2_2, mid_ijk_r2r_0_2, hstride);
          intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2_2, mid_ijk_r2r_0_2, mid_ijk_r2r_1_2, hstride);
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
