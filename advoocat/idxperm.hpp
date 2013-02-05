/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/
/** @mainpage
 *  TODO: table of examples with columns:
 *  - adv. algorithm
 *  - Jacobian
 *  - number of dimensions
 *  - equation set
 */

// code licensed under the terms of GNU GPL v3
// copyright holder: University of Warsaw

#pragma once

namespace advoocat
{
  namespace idxperm
  {
    // 2D

    template<int d> 
    inline idx_t<2> pi(
      const rng_t &i, 
      const rng_t &j
    );

    template<>
    inline idx_t<2> pi<0>(
      const rng_t &i, 
      const rng_t &j
    ) {
      return idx_t<2>({i,j});
    }

    template<>
    inline idx_t<2> pi<1>(
      const rng_t &j, 
      const rng_t &i
    ) {
      return idx_t<2>({i,j});
    } 

    // 3D

    template<int d> 
    inline idx_t<3> pi(
      const rng_t &i, 
      const rng_t &j,
      const rng_t &k
    );

    template<>
    inline idx_t<3> pi<0>(
      const rng_t &i, 
      const rng_t &j,
      const rng_t &k
    ) {
      return idx_t<3>({i,j,k});
    }

    template<>
    inline idx_t<3> pi<1>(
      const rng_t &j, 
      const rng_t &k,
      const rng_t &i
    ) {
      return idx_t<3>({i,j,k});
    }

    template<>
    inline idx_t<3> pi<2>(
      const rng_t &k, 
      const rng_t &i,
      const rng_t &j
    ) {
      return idx_t<3>({i,j,k});
    }
  }; // namespace idxperm
}; // namespace advoocat
