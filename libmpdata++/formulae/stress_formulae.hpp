/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/nabla_formulae.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace stress
    {
      // velocity gradient
      // 2D version
      template <int nd, class arrvec_t, class ijk_t, class dijk_t>
      inline void calc_vgrad(arrvec_t &vg,
                             const arrvec_t &v,
                             const ijk_t &ijk,
                             const dijk_t &dijk,
                             typename std::enable_if<nd == 2>::type* = 0)
      {
        vg[0](ijk) = formulae::nabla::grad<0>(v[0], ijk[0], ijk[1], dijk[0]);
        vg[1](ijk) = formulae::nabla::grad<0>(v[1], ijk[0], ijk[1], dijk[0]);

        vg[2](ijk) = formulae::nabla::grad<1>(v[0], ijk[1], ijk[0], dijk[1]);
        vg[3](ijk) = formulae::nabla::grad<1>(v[1], ijk[1], ijk[0], dijk[1]);
      }

      // 3D version
      template <int nd, class arrvec_t, class ijk_t, class dijk_t>
      inline void calc_vgrad(arrvec_t &vg,
                             const arrvec_t &v,
                             const ijk_t &ijk,
                             const dijk_t &dijk,
                             typename std::enable_if<nd == 3>::type* = 0)
      {
        vg[0](ijk) = formulae::nabla::grad<0>(v[0], ijk[0], ijk[1], ijk[2], dijk[0]);
        vg[1](ijk) = formulae::nabla::grad<0>(v[1], ijk[0], ijk[1], ijk[2], dijk[0]);
        vg[2](ijk) = formulae::nabla::grad<0>(v[2], ijk[0], ijk[1], ijk[2], dijk[0]);

        vg[3](ijk) = formulae::nabla::grad<1>(v[0], ijk[1], ijk[2], ijk[0], dijk[1]);
        vg[4](ijk) = formulae::nabla::grad<1>(v[1], ijk[1], ijk[2], ijk[0], dijk[1]);
        vg[5](ijk) = formulae::nabla::grad<1>(v[2], ijk[1], ijk[2], ijk[0], dijk[1]);

        vg[6](ijk) = formulae::nabla::grad<2>(v[0], ijk[2], ijk[0], ijk[1], dijk[2]);
        vg[7](ijk) = formulae::nabla::grad<2>(v[1], ijk[2], ijk[0], ijk[1], dijk[2]);
        vg[8](ijk) = formulae::nabla::grad<2>(v[2], ijk[2], ijk[0], ijk[1], dijk[2]);
      }
      
      // calculates unique deformation tensor components
      // 2D version
      template <int nd, class arrvec_t, class ijk_t>
      inline void calc_deform(arrvec_t &tau,
                              const arrvec_t &vg, 
                              const ijk_t &ijk,
                              typename std::enable_if<nd == 2>::type* = 0)
      {
        tau[0](ijk) = vg[0](ijk);
        tau[1](ijk) = 0.5 * (vg[1](ijk) + vg[2](ijk));
        tau[2](ijk) = vg[3](ijk);
      }
      
      // 3D version
      template <int nd, class arrvec_t, class ijk_t>
      inline void calc_deform(arrvec_t &tau,
                              const arrvec_t &vg,
                              const ijk_t &ijk,
                              typename std::enable_if<nd == 3>::type* = 0)
      {
        tau[0](ijk) = vg[0](ijk);
        tau[1](ijk) = 0.5 * (vg[3](ijk) + vg[1](ijk));
        tau[2](ijk) = 0.5 * (vg[6](ijk) + vg[2](ijk));
        tau[3](ijk) = vg[4](ijk);
        tau[4](ijk) = 0.5 * (vg[5](ijk) + vg[7](ijk));
        tau[5](ijk) = vg[8](ijk);
      }
      
      // calculate elements of stress tensor divergence
      // 2D version
      template <int nd, class arrvec_t, class ijk_t, class dijk_t>
      inline void calc_stress_div(arrvec_t &sdiv,
                                  const arrvec_t &tau,
                                  const ijk_t &ijk,
                                  const dijk_t &dijk,
                                  typename std::enable_if<nd == 2>::type* = 0)
      {
        sdiv[0](ijk) = formulae::nabla::grad<0>(tau[0], ijk[0], ijk[1], dijk[0]);
        sdiv[1](ijk) = formulae::nabla::grad<0>(tau[1], ijk[0], ijk[1], dijk[0]);

        sdiv[2](ijk) = formulae::nabla::grad<1>(tau[1], ijk[1], ijk[0], dijk[1]);
        sdiv[3](ijk) = formulae::nabla::grad<1>(tau[2], ijk[1], ijk[0], dijk[1]);
      }
      
      // 3D version
      template <int nd, class arrvec_t, class ijk_t, class dijk_t>
      inline void calc_stress_div(arrvec_t &sdiv,
                                  const arrvec_t &tau,
                                  const ijk_t &ijk,
                                  const dijk_t &dijk,
                                  typename std::enable_if<nd == 3>::type* = 0)
      {
        sdiv[0](ijk) = formulae::nabla::grad<0>(tau[0], ijk[0], ijk[1], ijk[2], dijk[0]);
        sdiv[1](ijk) = formulae::nabla::grad<0>(tau[1], ijk[0], ijk[1], ijk[2], dijk[0]);
        sdiv[2](ijk) = formulae::nabla::grad<0>(tau[2], ijk[0], ijk[1], ijk[2], dijk[0]);

        sdiv[3](ijk) = formulae::nabla::grad<1>(tau[1], ijk[1], ijk[2], ijk[0], dijk[1]);
        sdiv[4](ijk) = formulae::nabla::grad<1>(tau[3], ijk[1], ijk[2], ijk[0], dijk[1]);
        sdiv[5](ijk) = formulae::nabla::grad<1>(tau[4], ijk[1], ijk[2], ijk[0], dijk[1]);

        sdiv[6](ijk) = formulae::nabla::grad<2>(tau[2], ijk[2], ijk[0], ijk[1], dijk[2]);
        sdiv[7](ijk) = formulae::nabla::grad<2>(tau[4], ijk[2], ijk[0], ijk[1], dijk[2]);
        sdiv[8](ijk) = formulae::nabla::grad<2>(tau[5], ijk[2], ijk[0], ijk[1], dijk[2]);
      }

      // Pade correction
      using idxperm::pi;
      // 2D version
      template <int d, class arg_t>
      inline auto pade_helper(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j
      ) return_macro(,
	  x(idxperm::pi<d>(i + 1, j)) + 2 * x(idxperm::pi<d>(i, j)) + x(idxperm::pi<d>(i - 1, j))
      )

      template <int nd, class arrvec_t, class ijk_t>
      inline void pade_dispatch(arrvec_t &wrk, const ijk_t &ijk, int d, typename std::enable_if<nd == 2>::type* = 0)
      {
        switch(d)
        {
          case 0 : case 1 :
          {
            wrk[1](ijk) = pade_helper<0>(wrk[0], ijk[0], ijk[1]);
            break;
          }
          case 2 : case 3 :
          {
            wrk[1](ijk) = pade_helper<1>(wrk[0], ijk[1], ijk[0]);
            break;
          }
          default : assert(false);
        }
      }

      // 3D version
      template <int d, class arg_t>
      inline auto pade_helper(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k
      ) return_macro(,
	  x(idxperm::pi<d>(i+1, j, k)) + 2 * x(idxperm::pi<d>(i, j, k)) + x(idxperm::pi<d>(i - 1, j, k))
      )

      template <int nd, class arrvec_t, class ijk_t>
      inline void pade_dispatch(arrvec_t &wrk, const ijk_t &ijk, int d, typename std::enable_if<nd == 3>::type* = 0)
      {
        switch(d)
        {
          case 0 : case 1 : case 2 :
          {
            wrk[1](ijk) = pade_helper<0>(wrk[0], ijk[0], ijk[1], ijk[2]);
            break;
          }
          case 3 : case 4 : case 5 :
          {
            wrk[1](ijk) = pade_helper<1>(wrk[0], ijk[1], ijk[2], ijk[0]);
            break;
          }
          case 6 : case 7 : case 8 :
          {
            wrk[1](ijk) = pade_helper<2>(wrk[0], ijk[2], ijk[0], ijk[1]);
            break;
          }
          default : assert(false);
        }
      }

      // add stress forces
      // 2D version
      template <int nd, class arrvec_t, class ijk_t, class real_t>
      inline void calc_stress_rhs(arrvec_t &rhs, const arrvec_t &drv, ijk_t &ijk, real_t coeff, typename std::enable_if<nd == 2>::type* = 0)
      {
        rhs[0](ijk) += coeff * (drv[0](ijk) + drv[2](ijk));
        rhs[1](ijk) += coeff * (drv[1](ijk) + drv[3](ijk));
      }
      
      // 3D version
      template <int nd, class arrvec_t, class ijk_t, class real_t>
      inline void calc_stress_rhs(arrvec_t &rhs, const arrvec_t &drv, ijk_t &ijk, real_t coeff, typename std::enable_if<nd == 3>::type* = 0)
      {
        rhs[0](ijk) += coeff * (drv[0](ijk) + drv[3](ijk) + drv[6](ijk));
        rhs[1](ijk) += coeff * (drv[1](ijk) + drv[4](ijk) + drv[7](ijk));
        rhs[2](ijk) += coeff * (drv[2](ijk) + drv[5](ijk) + drv[8](ijk));
      }
    } // namespace stress
  } // namespace formulae
} // namespace libmpdataxx
