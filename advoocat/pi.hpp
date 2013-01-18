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

template<int d> 
inline idx_2d_t pi(const rng_t &i, const rng_t &j);

template<>
inline idx_2d_t pi<0>(
  const rng_t &i, const rng_t &j
) {
  return idx_2d_t({i,j});
};

template<>
inline idx_2d_t pi<1>(
  const rng_t &j, const rng_t &i
) {
  return idx_2d_t({i,j});
}; 
