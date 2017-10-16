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
      using arakawa_c::h;
      using opts::opts_t;

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
        tau[0](ijk) = 2 * vg[0](ijk);
        tau[1](ijk) = vg[1](ijk) + vg[2](ijk);
        tau[2](ijk) = 2 * vg[3](ijk);
      }
      
      // 3D version
      template <int nd, class arrvec_t, class ijk_t>
      inline void calc_deform(arrvec_t &tau,
                              const arrvec_t &vg,
                              const ijk_t &ijk,
                              typename std::enable_if<nd == 3>::type* = 0)
      {
        tau[0](ijk) = 2 * vg[0](ijk);
        tau[1](ijk) = vg[3](ijk) + vg[1](ijk);
        tau[2](ijk) = vg[6](ijk) + vg[2](ijk);
        tau[3](ijk) = 2 * vg[4](ijk);
        tau[4](ijk) = vg[5](ijk) + vg[7](ijk);
        tau[5](ijk) = 2 * vg[8](ijk);
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
      
      // Total deformation
      // 2D version
      template <int nd, class arrvec_t, class ijk_t>
      inline auto calc_tdef_sq(
	const arrvec_t &tau,
        const ijk_t &ijk,
        typename std::enable_if<nd == 2>::type* = 0
      ) return_macro(,
	  pow2(tau[0](ijk)) / 2 + pow2(tau[1](ijk)) + pow2(tau[2](ijk)) / 2
      )
      
      // 3D version
      template <int nd, class arrvec_t, class ijk_t>
      inline auto calc_tdef_sq(
	const arrvec_t &tau,
        const ijk_t &ijk,
        typename std::enable_if<nd == 3>::type* = 0
      ) return_macro(,
          pow2(tau[0](ijk)) / 2 + 
          pow2(tau[1](ijk)) +
          pow2(tau[2](ijk)) +
          pow2(tau[3](ijk)) / 2 +
          pow2(tau[4](ijk)) +
          pow2(tau[5](ijk)) / 2
      )
      
      // Compact formulation
      
      // surface stress
      // 2D version
      template <int nd, opts_t opts, class arr_t, class arrvec_t, class real_t, class ijk_t, class ijkm_t>
      inline void calc_drag_cmpct(arrvec_t &tau,
                                  const arrvec_t &v,
                                  const arr_t &rho,
                                  const real_t cdrag,
                                  const ijk_t &ijk,
                                  const ijkm_t &ijkm,
                                  typename std::enable_if<nd == 2>::type* = 0)
      {
        auto zro = rng_t(0, 0);
        tau[0](ijkm[0] + h, zro) = cdrag / 8 * 
                                   abs((v[0](ijkm[0] + 1, zro) + v[0](ijkm[0], zro))) *
                                   (v[0](ijkm[0] + 1, zro) + v[0](ijkm[0], zro)) *
                                   (  G<opts, 0>(rho, ijkm[0] + 1, zro)
                                    + G<opts, 0>(rho, ijkm[0]    , zro) );

      }

      // 3D version
      template <int nd, opts_t opts, class arr_t, class arrvec_t, class real_t, class ijk_t, class ijkm_t>
      inline void calc_drag_cmpct(arrvec_t &tau,
                                  const arrvec_t &v,
                                  const arr_t &rho,
                                  const real_t cdrag,
                                  const ijk_t &ijk,
                                  const ijkm_t &ijkm,
                                  typename std::enable_if<nd == 3>::type* = 0)
      {
        auto zro = rng_t(0, 0);
        tau[0](ijkm[0] + h, ijk[1], zro) = cdrag / 8 * sqrt(
                                                pow2((v[0](ijkm[0] + 1, ijk[1], zro) + v[0](ijkm[0], ijk[1], zro)))
                                              + pow2((v[1](ijkm[0] + 1, ijk[1], zro) + v[1](ijkm[0], ijk[1], zro)))
                                              ) *
                                           (v[0](ijkm[0] + 1, ijk[1], zro) + v[0](ijkm[0], ijk[1], zro)) *
                                           (  G<opts, 0>(rho, ijkm[0] + 1, ijk[1], zro)
                                            + G<opts, 0>(rho, ijkm[0]    , ijk[1], zro) );

        tau[1](ijk[0], ijkm[1] + h, zro) = cdrag / 8 * sqrt(
                                                pow2((v[0](ijk[0], ijkm[1] + 1, zro) + v[0](ijk[0], ijkm[1], zro)))
                                              + pow2((v[1](ijk[0], ijkm[1] + 1, zro) + v[1](ijk[0], ijkm[1], zro)))
                                              ) *
                                           (v[1](ijk[0], ijkm[1] + 1, zro) + v[1](ijk[0], ijkm[1], zro)) *
                                           (  G<opts, 0>(rho, ijk[0], ijkm[1] + 1, zro)
                                            + G<opts, 0>(rho, ijk[0], ijkm[1]    , zro) );
      }
      
      // velocity divergence
      // 2D version
      template <int nd, class arr_t, class arrvec_t, class ijk_t, class dijk_t>
      inline void calc_vip_div_cmpct(arr_t &div,
                                  const arrvec_t &v,
                                  const arr_t &rho,
                                  const ijk_t &ijk,
                                  const dijk_t &dijk,
                                  typename std::enable_if<nd == 2>::type* = 0)
      {
        div(ijk[0], ijk[1] + h) = - (v[1](ijk[0], ijk[1] + 1) + v[1](ijk[0], ijk[1])) *
                                    (rho(ijk[0], ijk[1] + 1) - rho(ijk[0], ijk[1])) /
                                    (dijk[1] * (rho(ijk[0], ijk[1] + 1) + rho(ijk[0], ijk[1])));
      }
      
      // 3D version
      template <int nd, class arr_t, class arrvec_t, class ijk_t, class dijk_t>
      inline void calc_vip_div_cmpct(arr_t &div,
                                  const arrvec_t &v,
                                  const arr_t &rho,
                                  const ijk_t &ijk,
                                  const dijk_t &dijk,
                                  typename std::enable_if<nd == 3>::type* = 0)
      {
        div(ijk[0], ijk[1], ijk[2] + h) = - (v[2](ijk[0], ijk[1], ijk[2] + 1) + v[2](ijk[0], ijk[1], ijk[2])) *
                                            (rho(ijk[0], ijk[1], ijk[2] + 1) - rho(ijk[0], ijk[1], ijk[2])) /
                                            (dijk[2] * (rho(ijk[0], ijk[1], ijk[2] + 1) + rho(ijk[0], ijk[1], ijk[2])));
      }

      // calculates unique deformation tensor components
      // 2D version
      template <int nd, class arrvec_t, class arr_t, class ijk_t, class ijkm_t, class dijk_t>
      inline void calc_deform_cmpct(arrvec_t &tau,
                                    const arrvec_t &v,
                                    const arr_t &div_v,
                                    const ijk_t &ijk,
                                    const ijkm_t &ijkm,
                                    const dijk_t &dijk,
                                    typename std::enable_if<nd == 2>::type* = 0)
      {
        tau[0](ijkm[0] + h, ijk[1])
        =
        2 * (
          (v[0](ijkm[0] + 1, ijk[1]) - v[0](ijkm[0], ijk[1])) / dijk[0]
          - 0.25 / 2 * ( div_v(ijkm[0], ijk[1] + h) + div_v(ijkm[0], ijk[1] - h) 
                       + div_v(ijkm[0] + 1, ijk[1] + h) + div_v(ijkm[0] + 1, ijk[1] - h))
        );
        
        tau[1](ijk[0], ijkm[1] + h)
        =
        2 * (
          (v[1](ijk[0], ijkm[1] + 1) - v[1](ijk[0], ijkm[1])) / dijk[1]
          - 1.0 / 2 * div_v(ijk[0], ijkm[1] + h)
        );

        tau[2](ijkm[0] + h, ijkm[1] + h)
        =
        0.5 *
        (
          ( v[0](ijkm[0], ijkm[1] + 1) - v[0](ijkm[0], ijkm[1]) 
          + v[0](ijkm[0] + 1, ijkm[1] + 1) - v[0](ijkm[0] + 1, ijkm[1])
          ) / dijk[1]
          +
          ( v[1](ijkm[0] + 1, ijkm[1]) - v[1](ijkm[0], ijkm[1]) 
          + v[1](ijkm[0] + 1, ijkm[1] + 1) - v[1](ijkm[0], ijkm[1] + 1)
          ) / dijk[0]
        );
      }

      // 3D version
      template <int nd, class arrvec_t, class arr_t, class ijk_t, class ijkm_t, class dijk_t>
      inline void calc_deform_cmpct(arrvec_t &tau,
                                    const arrvec_t &v,
                                    const arr_t &div_v,
                                    const ijk_t &ijk,
                                    const ijkm_t &ijkm,
                                    const dijk_t &dijk,
                                    typename std::enable_if<nd == 3>::type* = 0)
      {
        tau[0](ijkm[0] + h, ijk[1], ijk[2])
        =
        2 * (
          (v[0](ijkm[0] + 1, ijk[1], ijk[2]) - v[0](ijkm[0], ijk[1], ijk[2])) / dijk[0]
          - 0.25 / 3 * ( div_v(ijkm[0], ijk[1], ijk[2] + h) + div_v(ijkm[0], ijk[1], ijk[2] - h)
                       + div_v(ijkm[0] + 1, ijk[1], ijk[2] + h) + div_v(ijkm[0] + 1, ijk[1], ijk[2] - h))
        );
        
        tau[1](ijk[0], ijkm[1] + h, ijk[2])
        =
        2 * (
          (v[1](ijk[0], ijkm[1] + 1, ijk[2]) - v[1](ijk[0], ijkm[1], ijk[2])) / dijk[1]
          - 0.25 / 3 * ( div_v(ijk[0], ijkm[1], ijk[2] + h) + div_v(ijk[0], ijkm[1], ijk[2] - h)
                       + div_v(ijk[0], ijkm[1] + 1, ijk[2] + h) + div_v(ijk[0], ijkm[1] + 1, ijk[2] - h))
        );
        
        tau[2](ijk[0], ijk[1], ijkm[2] + h)
        =
        2 * (
          (v[2](ijk[0], ijk[1], ijkm[2] + 1) - v[2](ijk[0], ijk[1], ijkm[2])) / dijk[2]
          - 1.0 / 3  * div_v(ijk[0], ijk[1], ijkm[2] + h) 
        );

        tau[3](ijkm[0] + h, ijkm[1] + h, ijk[2])
        =
        0.5 *
        (
          ( v[0](ijkm[0], ijkm[1] + 1, ijk[2]) - v[0](ijkm[0], ijkm[1], ijk[2]) 
          + v[0](ijkm[0] + 1, ijkm[1] + 1, ijk[2]) - v[0](ijkm[0] + 1, ijkm[1], ijk[2])
          ) / dijk[1]
          +
          ( v[1](ijkm[0] + 1, ijkm[1], ijk[2]) - v[1](ijkm[0], ijkm[1], ijk[2]) 
          + v[1](ijkm[0] + 1, ijkm[1] + 1, ijk[2]) - v[1](ijkm[0], ijkm[1] + 1, ijk[2])
          ) / dijk[0]
        );
        
        tau[4](ijkm[0] + h, ijk[1], ijkm[2] + h)
        =
        0.5 *
        (
          ( v[0](ijkm[0], ijk[1], ijkm[2] + 1) - v[0](ijkm[0], ijk[1], ijkm[2])
          + v[0](ijkm[0] + 1, ijk[1], ijkm[2] + 1) - v[0](ijkm[0] + 1, ijk[1], ijkm[2])
          ) / dijk[2]
          +
          ( v[2](ijkm[0] + 1, ijk[1], ijkm[2]) - v[2](ijkm[0], ijk[1], ijkm[2]) 
          + v[2](ijkm[0] + 1, ijk[1], ijkm[2] + 1) - v[2](ijkm[0], ijk[1], ijkm[2] + 1)
          ) / dijk[0]
        );
        
        tau[5](ijk[0], ijkm[1] + h, ijkm[2] + h)
        =
        0.5 *
        (
          ( v[1](ijk[0], ijkm[1], ijkm[2] + 1) - v[1](ijk[0], ijkm[1], ijkm[2])
          + v[1](ijk[0], ijkm[1] + 1, ijkm[2] + 1) - v[1](ijk[0], ijkm[1] + 1, ijkm[2])
          ) / dijk[2]
          +
          ( v[2](ijk[0], ijkm[1] + 1, ijkm[2]) - v[2](ijk[0], ijkm[1], ijkm[2]) 
          + v[2](ijk[0], ijkm[1] + 1, ijkm[2] + 1) - v[2](ijk[0], ijkm[1], ijkm[2] + 1)
          ) / dijk[1]
        );
      }
      
      // Total deformation
      // 2D version
      template <int nd, class arrvec_t, class ijk_t>
      inline auto calc_tdef_sq_cmpct(
	const arrvec_t &tau,
        const ijk_t &ijk,
        typename std::enable_if<nd == 2>::type* = 0
      ) return_macro(
        ,
          // one half taken as an average 
          (
            pow2(tau[0](ijk[0] + h, ijk[1])) + pow2(tau[0](ijk[0] - h, ijk[1]))
          + pow2(tau[1](ijk[0], ijk[1] + h)) + pow2(tau[1](ijk[0], ijk[1] - h))
          + pow2(tau[2](ijk[0] + h, ijk[1] + h)) + pow2(tau[2](ijk[0] - h, ijk[1] + h))
          + pow2(tau[2](ijk[0] + h, ijk[1] - h)) + pow2(tau[2](ijk[0] - h, ijk[1] - h))
          ) / 8
          +
          // second half taken as max TODO: as an option ?
          (
            max<rng_t>(pow2(tau[0](ijk[0] + h, ijk[1])), pow2(tau[0](ijk[0] - h, ijk[1])))
          + max<rng_t>(pow2(tau[1](ijk[0], ijk[1] + h)), pow2(tau[1](ijk[0], ijk[1] - h)))
          + 2 * max<rng_t>(pow2(tau[2](ijk[0] + h, ijk[1] + h)), pow2(tau[2](ijk[0] - h, ijk[1] + h)),
                           pow2(tau[2](ijk[0] + h, ijk[1] - h)), pow2(tau[2](ijk[0] - h, ijk[1] - h)))
          ) / 4
      )

      // 3D version
      template <int nd, class arrvec_t, class ijk_t>
      inline auto calc_tdef_sq_cmpct(
	const arrvec_t &tau,
        const ijk_t &ijk,
        typename std::enable_if<nd == 3>::type* = 0
      ) return_macro(
        ,
          // one half taken as an average 
          (
            pow2(tau[0](ijk[0] + h, ijk[1], ijk[2])) + pow2(tau[0](ijk[0] - h, ijk[1], ijk[2]))
          + pow2(tau[1](ijk[0], ijk[1] + h, ijk[2])) + pow2(tau[1](ijk[0], ijk[1] - h, ijk[2]))
          + pow2(tau[2](ijk[0], ijk[1], ijk[2] + h)) + pow2(tau[2](ijk[0], ijk[1], ijk[2] - h))
          + pow2(tau[3](ijk[0] + h, ijk[1] + h, ijk[2])) + pow2(tau[3](ijk[0] - h, ijk[1] + h, ijk[2]))
          + pow2(tau[3](ijk[0] + h, ijk[1] - h, ijk[2])) + pow2(tau[3](ijk[0] - h, ijk[1] - h, ijk[2]))
          + pow2(tau[4](ijk[0] + h, ijk[1], ijk[2] + h)) + pow2(tau[4](ijk[0] - h, ijk[1], ijk[2] + h))
          + pow2(tau[4](ijk[0] + h, ijk[1], ijk[2] - h)) + pow2(tau[4](ijk[0] - h, ijk[1], ijk[2] - h))
          + pow2(tau[5](ijk[0], ijk[1] + h, ijk[2] + h)) + pow2(tau[5](ijk[0], ijk[1] - h, ijk[2] + h))
          + pow2(tau[5](ijk[0], ijk[1] + h, ijk[2] - h)) + pow2(tau[5](ijk[0], ijk[1] - h, ijk[2] - h))
          ) / 8
          +
          // second half taken as max TODO: as an option ?
          (
            max<rng_t>(pow2(tau[0](ijk[0] + h, ijk[1], ijk[2])), pow2(tau[0](ijk[0] - h, ijk[1], ijk[2])))
          + max<rng_t>(pow2(tau[1](ijk[0], ijk[1] + h, ijk[2])), pow2(tau[1](ijk[0], ijk[1] - h, ijk[2])))
          + max<rng_t>(pow2(tau[2](ijk[0], ijk[1], ijk[2] + h)), pow2(tau[2](ijk[0], ijk[1], ijk[2] - h)))
          + 2 * max<rng_t>(pow2(tau[3](ijk[0] + h, ijk[1] + h, ijk[2])), pow2(tau[3](ijk[0] - h, ijk[1] + h, ijk[2])),
                           pow2(tau[3](ijk[0] + h, ijk[1] - h, ijk[2])), pow2(tau[3](ijk[0] - h, ijk[1] - h, ijk[2])))
          + 2 * max<rng_t>(pow2(tau[4](ijk[0] + h, ijk[1], ijk[2] + h)), pow2(tau[4](ijk[0] - h, ijk[1], ijk[2] + h)),
                           pow2(tau[4](ijk[0] + h, ijk[1], ijk[2] - h)), pow2(tau[4](ijk[0] - h, ijk[1], ijk[2] - h)))
          + 2 * max<rng_t>(pow2(tau[5](ijk[0], ijk[1] + h, ijk[2] + h)), pow2(tau[5](ijk[0], ijk[1] - h, ijk[2] + h)),
                           pow2(tau[5](ijk[0], ijk[1] + h, ijk[2] - h)), pow2(tau[5](ijk[0], ijk[1] - h, ijk[2] - h)))
          ) / 4
      )
     
      // multiplication of compact vector components by constant molecular viscosity
      // 2D version
      template <int nd, class arrvec_t, class real_t, class ijk_t>
      inline void multiply_vctr_cmpct(const arrvec_t &av,
                                     real_t coeff,
                                     const ijk_t &ijk,
                                     typename std::enable_if<nd == 2>::type* = 0)
      {
        av[0](ijk[0] + h, ijk[1]) *= coeff;
        av[1](ijk[0], ijk[1] + h) *= coeff;
      }

      // 3D version
      template <int nd, class arrvec_t, class real_t, class ijk_t>
      inline void multiply_vctr_cmpct(const arrvec_t &av,
                                     real_t coeff,
                                     const ijk_t &ijk,
                                     typename std::enable_if<nd == 3>::type* = 0)
      {
        av[0](ijk[0] + h, ijk[1], ijk[2]) *= coeff;
        av[1](ijk[0], ijk[1] + h, ijk[2]) *= coeff;
        av[2](ijk[0], ijk[1], ijk[2] + h) *= coeff;
      }

      // multiplication of compact vector components by variable eddy viscosity
      // 2D version
      template <int nd, opts_t opts, class arrvec_t, class real_t, class arr_t, class ijk_t>
      inline void multiply_vctr_cmpct(const arrvec_t &av,
                                      real_t coeff,
                                      const arr_t & k_m,
                                      const arr_t & rho,
                                      const ijk_t &ijk,
                                      typename std::enable_if<nd == 2>::type* = 0)
      {
        av[0](ijk[0] + h, ijk[1]) *= coeff * 
                           0.5 * (k_m(ijk[0] + 1, ijk[1]) + k_m(ijk[0], ijk[1])) *
                           0.5 * (G<opts, 0>(rho, ijk[0] + 1, ijk[1]) + G<opts, 0>(rho, ijk[0], ijk[1]));
        
        av[1](ijk[0], ijk[1] + h) *= coeff *
                           0.5 * (k_m(ijk[0], ijk[1] + 1) + k_m(ijk[0], ijk[1])) * 
                           0.5 * (G<opts, 0>(rho, ijk[0], ijk[1] + 1) + G<opts, 0>(rho, ijk[0], ijk[1]));
      }

      // 3D version
      template <int nd, opts_t opts, class arrvec_t, class real_t, class arr_t, class ijk_t>
      inline void multiply_vctr_cmpct(const arrvec_t &av,
                                      real_t coeff,
                                      const arr_t & k_m,
                                      const arr_t & rho,
                                      const ijk_t &ijk,
                                      typename std::enable_if<nd == 3>::type* = 0)
      {
        av[0](ijk[0] + h, ijk[1], ijk[2]) *= coeff *
                           0.5 * (k_m(ijk[0] + 1, ijk[1], ijk[2]) + k_m(ijk[0], ijk[1], ijk[2])) *
                           0.5 * (G<opts, 0>(rho, ijk[0] + 1, ijk[1], ijk[2]) + G<opts, 0>(rho,ijk[0], ijk[1], ijk[2]));
        
        av[1](ijk[0], ijk[1] + h, ijk[2]) *= coeff *
                           0.5 * (k_m(ijk[0], ijk[1] + 1, ijk[2]) + k_m(ijk[0], ijk[1], ijk[2])) *
                           0.5 * (G<opts, 0>(rho, ijk[0], ijk[1] + 1, ijk[2]) + G<opts, 0>(rho, ijk[0], ijk[1], ijk[2]));
        
        av[2](ijk[0], ijk[1], ijk[2] + h) *= coeff *
                           0.5 * (k_m(ijk[0], ijk[1], ijk[2] + 1) + k_m(ijk[0], ijk[1], ijk[2])) *
                           0.5 * (G<opts, 0>(rho, ijk[0], ijk[1], ijk[2] + 1) + G<opts, 0>(rho, ijk[0], ijk[1], ijk[2]));
      }
      
      // multiplication of compact tensor components by constant molecular viscosity
      // 2D version
      template <int nd, class arrvec_t, class real_t, class ijk_t>
      inline void multiply_tnsr_cmpct(const arrvec_t &av,
                                      const real_t coeff,
                                      const ijk_t &ijk,
                                      typename std::enable_if<nd == 2>::type* = 0)
      {
        multiply_vctr_cmpct<nd>(av, coeff, ijk);
        av[2](ijk[0] + h, ijk[1] + h) *= coeff;
      }

      template <int nd, class arrvec_t, class real_t, class ijk_t>
      inline void multiply_tnsr_cmpct(const arrvec_t &av,
                                      const real_t coeff,
                                      const ijk_t &ijk,
                                      typename std::enable_if<nd == 3>::type* = 0)
      {
        multiply_vctr_cmpct<nd>(av, coeff, ijk);
        av[3](ijk[0] + h, ijk[1] + h, ijk[2]) *= coeff;
        av[4](ijk[0] + h, ijk[1], ijk[2] + h) *= coeff;
        av[5](ijk[0], ijk[1] + h, ijk[2] + h) *= coeff;
      }

      // multiplication of compact tensor components by variable eddy viscosity
      // 2D version
      template <int nd, opts_t opts, class arrvec_t, class real_t, class arr_t, class ijk_t>
      inline void multiply_tnsr_cmpct(const arrvec_t &av,
                                      const real_t coeff,
                                      const arr_t &k_m,
                                      const arr_t &rho,
                                      const ijk_t &ijk,
                                      typename std::enable_if<nd == 2>::type* = 0)
      {
        multiply_vctr_cmpct<nd, opts>(av, coeff, k_m, rho, ijk);
        av[2](ijk[0] + h, ijk[1] + h) *= coeff * 
                                         0.25 * ( k_m(ijk[0] + 1, ijk[1]    )
                                                + k_m(ijk[0]    , ijk[1]    )
                                                + k_m(ijk[0] + 1, ijk[1] + 1)
                                                + k_m(ijk[0]    , ijk[1] + 1)
                                                ) *
                                         0.25 * ( G<opts, 0>(rho, ijk[0] + 1, ijk[1]    )
                                                + G<opts, 0>(rho, ijk[0]    , ijk[1]    )
                                                + G<opts, 0>(rho, ijk[0] + 1, ijk[1] + 1)
                                                + G<opts, 0>(rho, ijk[0]    , ijk[1] + 1)
                                                );
      }

      // 3D version
      template <int nd, opts_t opts, class arrvec_t, class real_t, class arr_t, class ijk_t>
      inline void multiply_tnsr_cmpct(const arrvec_t &av,
                                      const real_t coeff,
                                      const arr_t &k_m,
                                      const arr_t &rho,
                                      const ijk_t &ijk,
                                      typename std::enable_if<nd == 3>::type* = 0)
      {
        multiply_vctr_cmpct<nd, opts>(av, coeff, k_m, rho, ijk);
        av[3](ijk[0] + h, ijk[1] + h, ijk[2]) *= coeff *
                                                 0.25 * ( k_m(ijk[0] + 1, ijk[1]    , ijk[2])
                                                        + k_m(ijk[0]    , ijk[1]    , ijk[2])
                                                        + k_m(ijk[0] + 1, ijk[1] + 1, ijk[2])
                                                        + k_m(ijk[0]    , ijk[1] + 1, ijk[2])
                                                        ) *
                                                 0.25 * ( G<opts, 0>(rho, ijk[0] + 1, ijk[1]    , ijk[2])
                                                        + G<opts, 0>(rho, ijk[0]    , ijk[1]    , ijk[2])
                                                        + G<opts, 0>(rho, ijk[0] + 1, ijk[1] + 1, ijk[2])
                                                        + G<opts, 0>(rho, ijk[0]    , ijk[1] + 1, ijk[2])
                                                        );
        
        av[4](ijk[0] + h, ijk[1], ijk[2] + h) *= coeff *
                                                 0.25 * ( k_m(ijk[0] + 1, ijk[1], ijk[2]    )
                                                        + k_m(ijk[0]    , ijk[1], ijk[2]    )
                                                        + k_m(ijk[0] + 1, ijk[1], ijk[2] + 1)
                                                        + k_m(ijk[0]    , ijk[1], ijk[2] + 1)
                                                        ) *
                                                 0.25 * ( G<opts, 0>(rho, ijk[0] + 1, ijk[1], ijk[2]    )
                                                        + G<opts, 0>(rho, ijk[0]    , ijk[1], ijk[2]    )
                                                        + G<opts, 0>(rho, ijk[0] + 1, ijk[1], ijk[2] + 1)
                                                        + G<opts, 0>(rho, ijk[0]    , ijk[1], ijk[2] + 1)
                                                        );

        av[5](ijk[0], ijk[1] + h, ijk[2] + h) *= coeff *
                                                 0.25 * ( k_m(ijk[0], ijk[1]    , ijk[2] + 1)
                                                        + k_m(ijk[0], ijk[1]    , ijk[2]    )
                                                        + k_m(ijk[0], ijk[1] + 1, ijk[2] + 1)
                                                        + k_m(ijk[0], ijk[1] + 1, ijk[2]    )
                                                        ) *
                                                 0.25 * ( G<opts, 0>(rho, ijk[0], ijk[1]    , ijk[2] + 1)
                                                        + G<opts, 0>(rho, ijk[0], ijk[1]    , ijk[2]    )
                                                        + G<opts, 0>(rho, ijk[0], ijk[1] + 1, ijk[2] + 1)
                                                        + G<opts, 0>(rho, ijk[0], ijk[1] + 1, ijk[2]    )
                                                        );
      }

      // flux divergence
      // 2D version
      template <int nd, opts_t opts, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline auto flux_div_cmpct(
	const arrvec_t &f,
	const arr_t &rho,
	const ijk_t &ijk,
	const dijk_t dijk,
        typename std::enable_if<nd == 2>::type* = 0
      ) return_macro(,
	( (f[0](ijk[0]+h, ijk[1]) - f[0](ijk[0]-h, ijk[1])) / dijk[0]
	  +
	  (f[1](ijk[0], ijk[1]+h) - f[1](ijk[0], ijk[1]-h)) / dijk[1]
        ) / G<opts>(rho, ijk)
      )
      
      // 3D version
      template <int nd, opts_t opts, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline auto flux_div_cmpct(
	const arrvec_t &f,
	const arr_t &rho,
	const ijk_t &ijk,
	const dijk_t dijk,
        typename std::enable_if<nd == 3>::type* = 0
      ) return_macro(,
	( (f[0](ijk[0]+h, ijk[1], ijk[2]) - f[0](ijk[0]-h, ijk[1], ijk[2])) / dijk[0]
	  +
	  (f[1](ijk[0], ijk[1]+h, ijk[2]) - f[1](ijk[0], ijk[1]-h, ijk[2])) / dijk[1]
	  +
	  (f[2](ijk[0], ijk[1], ijk[2]+h) - f[2](ijk[0], ijk[1], ijk[2]-h)) / dijk[2]
        ) / G<opts>(rho, ijk)
      )
     
      // add stress forces
      // 2D version
      template <int nd, opts_t opts, class arrvec_t, class arr_t, class ijk_t, class dijk_t, class real_t>
      inline void calc_stress_rhs_cmpct(arrvec_t &rhs,
                                        const arrvec_t &tau,
                                        const arr_t &rho,
                                        const ijk_t &ijk,
                                        const dijk_t &dijk,
                                        real_t coeff,
                                        typename std::enable_if<nd == 2>::type* = 0)
      {
        rhs[0](ijk)
        +=
        coeff *
        (
          (tau[0](ijk[0] + h, ijk[1]) - tau[0](ijk[0] - h, ijk[1])) / dijk[0]
          + 0.5 * 
          ( tau[2](ijk[0] + h, ijk[1] + h) - tau[2](ijk[0] + h, ijk[1] - h)
          + tau[2](ijk[0] - h, ijk[1] + h) - tau[2](ijk[0] - h, ijk[1] - h)
          ) / dijk[1]
        ) / G<opts>(rho, ijk);
        
        rhs[1](ijk)
        +=
        coeff *
        (
          0.5 * 
          ( tau[2](ijk[0] + h, ijk[1] + h) - tau[2](ijk[0] - h, ijk[1] + h)
          + tau[2](ijk[0] + h, ijk[1] - h) - tau[2](ijk[0] - h, ijk[1] - h)
          ) / dijk[0]
          +
          (tau[1](ijk[0], ijk[1] + h) - tau[1](ijk[0], ijk[1] - h)) / dijk[1]
        ) / G<opts>(rho, ijk);
      }

      // 3D version
      template <int nd, opts_t opts, class arrvec_t, class arr_t, class ijk_t, class dijk_t, class real_t>
      inline void calc_stress_rhs_cmpct(arrvec_t &rhs,
                                        const arrvec_t &tau,
                                        const arr_t &rho,
                                        const ijk_t &ijk,
                                        const dijk_t &dijk,
                                        real_t coeff,
                                        typename std::enable_if<nd == 3>::type* = 0)
      {
        rhs[0](ijk)
        +=
        coeff *
        (
          (tau[0](ijk[0] + h, ijk[1], ijk[2]) - tau[0](ijk[0] - h, ijk[1], ijk[2])) / dijk[0]
          + 0.5 * 
          ( tau[3](ijk[0] + h, ijk[1] + h, ijk[2]) - tau[3](ijk[0] + h, ijk[1] - h, ijk[2])
          + tau[3](ijk[0] - h, ijk[1] + h, ijk[2]) - tau[3](ijk[0] - h, ijk[1] - h, ijk[2])
          ) / dijk[1]
          + 0.5 * 
          ( tau[4](ijk[0] + h, ijk[1], ijk[2] + h) - tau[4](ijk[0] + h, ijk[1], ijk[2] - h)
          + tau[4](ijk[0] - h, ijk[1], ijk[2] + h) - tau[4](ijk[0] - h, ijk[1], ijk[2] - h)
          ) / dijk[2]
        ) / G<opts>(rho, ijk);
        
        rhs[1](ijk)
        +=
        coeff *
        (
          0.5 * 
          ( tau[3](ijk[0] + h, ijk[1] + h, ijk[2]) - tau[3](ijk[0] - h, ijk[1] + h, ijk[2])
          + tau[3](ijk[0] + h, ijk[1] - h, ijk[2]) - tau[3](ijk[0] - h, ijk[1] - h, ijk[2])
          ) / dijk[0]
          +
          (tau[1](ijk[0], ijk[1] + h, ijk[2]) - tau[1](ijk[0], ijk[1] - h, ijk[2])) / dijk[1]
          + 0.5 * 
          ( tau[5](ijk[0], ijk[1] + h, ijk[2] + h) - tau[5](ijk[0], ijk[1] + h, ijk[2] - h)
          + tau[5](ijk[0], ijk[1] - h, ijk[2] + h) - tau[5](ijk[0], ijk[1] - h, ijk[2] - h)
          ) / dijk[2]
        ) / G<opts>(rho, ijk);
        
        rhs[2](ijk)
        +=
        coeff *
        (
          0.5 * 
          ( tau[4](ijk[0] + h, ijk[1], ijk[2] + h) - tau[4](ijk[0] - h, ijk[1], ijk[2] + h)
          + tau[4](ijk[0] + h, ijk[1], ijk[2] - h) - tau[4](ijk[0] - h, ijk[1], ijk[2] - h)
          ) / dijk[0]
          + 0.5 * 
          ( tau[5](ijk[0], ijk[1] + h, ijk[2] + h) - tau[5](ijk[0], ijk[1] - h, ijk[2] + h)
          + tau[5](ijk[0], ijk[1] + h, ijk[2] - h) - tau[5](ijk[0], ijk[1] - h, ijk[2] - h)
          ) / dijk[1]
          +
          (tau[2](ijk[0], ijk[1], ijk[2] + h) - tau[2](ijk[0], ijk[1], ijk[2] - h)) / dijk[2]
        ) / G<opts>(rho, ijk);
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

    } // namespace stress
  } // namespace formulae
} // namespace libmpdataxx
