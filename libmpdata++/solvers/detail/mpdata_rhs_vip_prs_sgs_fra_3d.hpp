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
          arr_t &dst,
          const arr_t &src,
          const rng_t &i,
          const rng_t &j,
          const rng_t &k,
          const int &dist
        )
        {
          using idxperm::pi;
          using namespace arakawa_c;
          using real_t = typename ct_params_t::real_t;

          dst(pi<d>(i, j, k)) = real_t(.5) * (
            src(pi<d>(i - dist, j, k)) +
            src(pi<d>(i + dist, j, k))
          );
        }

        void interpolate_refinee(arrvec_t<typename parent_t::arr_t> &dst,
                                 const arrvec_t<typename parent_t::arr_t> &src,
                                 const int stride)
        {
          using namespace arakawa_c; // for rng_t operator^
          assert(dist % 2 == 0);
          const auto hstride = stride / 2;

//          auto ex = this->halo - 1;
          intrp<0>(dst[0], src[0], this->ijk_r2r[0]^hstride, this->ijk_r2r[1], this->ijk_r2r[2], hstride);
          intrp<1>(dst[1], src[1], this->ijk_r2r[1]^hstride, this->ijk_r2r[2], this->ijk_r2r[0], hstride);
          intrp<2>(dst[2], src[2], this->ijk_r2r[2]^hstride, this->ijk_r2r[0], this->ijk_r2r[1], hstride);
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
