/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

namespace libmpdataxx
{
  namespace idxperm
  {
    template <int n_dims> using int_idx_t = blitz::TinyVector<int, n_dims>;

    // 2D - ranges or mix of ranges with ints
    template<int d>
    inline idx_t<2> pi(const rng_t &i, const rng_t &j);

    template<>
    inline idx_t<2> pi<0>(const rng_t &i, const rng_t &j) { return idx_t<2>({i,j}); }

    template<>
    inline idx_t<2> pi<1>(const rng_t &j, const rng_t &i) { return idx_t<2>({i,j}); }

    template<int d>
    inline idx_t<2> pi(const int   &i, const rng_t &j) { return pi<d>(rng_t(i,i), j); }

    // 2D - ints
    template<int d>
    inline int_idx_t<2> pi(const int i, const int j);

    template<>
    inline int_idx_t<2> pi<0>(const int i, const int j) { return int_idx_t<2>({i,j}); }

    template<>
    inline int_idx_t<2> pi<1>(const int j, const int i) { return int_idx_t<2>({i,j}); }

    // 3D - ranges or mix of ranges with ints
    template<int d>
    inline idx_t<3> pi(const rng_t &i, const rng_t &j, const rng_t &k);

    template<>
    inline idx_t<3> pi<0>(const rng_t &i, const rng_t &j, const rng_t &k) { return idx_t<3>({i,j,k}); }

    template<>
    inline idx_t<3> pi<1>(const rng_t &j, const rng_t &k, const rng_t &i) { return idx_t<3>({i,j,k}); }

    template<>
    inline idx_t<3> pi<2>(const rng_t &k, const rng_t &i, const rng_t &j) { return idx_t<3>({i,j,k}); }

    template<int d>
    inline idx_t<3> pi(const int &i, const rng_t &j, const rng_t &k) { return pi<d>(rng_t(i,i), j, k); }

    template<int d>
    inline idx_t<3> pi(const int &i, const rng_t &j, const int   &k) { return pi<d>(rng_t(i,i), j, rng_t(k,k)); }

    template<int d>
    inline idx_t<3> pi(const int &i, const int   &j, const rng_t &k) { return pi<d>(rng_t(i,i), rng_t(j,j), k); }

    // 3D - ints
    template<int d>
    inline int_idx_t<3> pi(const int i, const int j, const int k);

    template<>
    inline int_idx_t<3> pi<0>(const int i, const int j, const int k) { return int_idx_t<3>({i,j,k}); }

    template<>
    inline int_idx_t<3> pi<1>(const int j, const int k, const int i) { return int_idx_t<3>({i,j,k}); }

    template<>
    inline int_idx_t<3> pi<2>(const int k, const int i, const int j) { return int_idx_t<3>({i,j,k}); }

    // 3D - strided ranges
    template<int d>
    inline idxs_t<3> pis(const rng_t &i, const rng_t &j, const rng_t &k);

    template<>
    inline idxs_t<3> pis<0>(const rng_t &i, const rng_t &j, const rng_t &k) { return idxs_t<3>{{i.first(), j.first(), k.first()}, {i.last(), j.last(), k.last()}, {i.stride(), j.stride(), k.stride()}}; }

    template<>
    inline idxs_t<3> pis<1>(const rng_t &j, const rng_t &k, const rng_t &i) { return idxs_t<3>{{i.first(), j.first(), k.first()}, {i.last(), j.last(), k.last()}, {i.stride(), j.stride(), k.stride()}}; }

    template<>
    inline idxs_t<3> pis<2>(const rng_t &k, const rng_t &i, const rng_t &j) { return idxs_t<3>{{i.first(), j.first(), k.first()}, {i.last(), j.last(), k.last()}, {i.stride(), j.stride(), k.stride()}}; }
  } // namespace idxperm
} // namespace libmpdataxx
