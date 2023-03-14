/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace fractal
    {
      // CDF of stretching parameter d; assuming CDF = A d^B + C under condition that CDF(0.5)=0 and CDF(1)=1
      template<class real_t>
      static real_t CDF_of_d(const real_t &d)
      {
        const real_t B = -0.39091161; // comes from a fit to the E. Akinlabi data
        return (pow(d,B) - pow(0.5,B)) / (1. - pow(0.5,B));
      }
      // inverse function: d of a CDF
      template<class real_t>
      static real_t d_of_CDF(const real_t &CDF)
      {
        const real_t B = -0.39091161; // comes from a fit to the E. Akinlabi data
        return pow((1.-pow(0.5,B))*CDF + pow(0.5,B), 1./B);
      }

      template<class real_t>
      struct d_of_CDF_fctr
      {
        real_t operator()(const real_t &CDF) const
        {
          return CDF < 0 ? -d_of_CDF(-CDF) : d_of_CDF(CDF);
        }
        BZ_DECLARE_FUNCTOR(d_of_CDF_fctr);
      };

      // (tri-)linear interpolation
      // 3d version
      template<int d, class real_t, class arr_t, class rng_t>
      void intrp(
        arr_t arr,       // array to be interpolated
        const rng_t &i,  // index of the point to be interpolated to
        const rng_t &j,  // ditto
        const rng_t &k,  // ditto
        const int &dist  // number of indices between an interpolated point and a known point
      )
      {
        using idxperm::pis;

        // debug output
         // std::cerr << "ranges in interpolation" << std::endl;
         // if(d==0)
         //   std::cerr << "range<" << d << ">: " << i << " " << j << " " << k << std::endl;
         // if(d==1)
         //   std::cerr << "range<" << d << ">: " << k << " " << i << " " << j << std::endl;
         // if(d==2)
         //   std::cerr << "range<" << d << ">: " << j << " " << k << " " << i << std::endl;

        arr(pis<d>(i, j, k)) = real_t(.5) * (
          arr(pis<d>(i - dist, j, k)) +
          arr(pis<d>(i + dist, j, k))
        );
      }

      // fractral reconstruction (calculates 2 mid-points based on 3)
      // following Akinlabi et al. "Fractal reconstruciton..." 2019
      // 3d version
      template<int d, class real_t, class arr_t, class rng_t>
      void rcnstrct(
        arr_t arr,       // array to be reconstructed
        const rng_t &i,  // index of the first point from the pair to be reconstructed
        const rng_t &j,  // ditto
        const rng_t &k,  // ditto
        const int &dist, // number of indices between a reconstructed point and a known point
        arr_t c_j,       // parameter c
        const arr_t d_j, // parameter d (stretching) - needs to be calculated before
        arr_t f_j        // parameter f
      )
      {
        // debug output
//          std::cerr << "ranges in rcnstrct" << std::endl;
//          if(d==0)
//            std::cerr << "range<" << d << ">: " << i << " " << j << " " << k << std::endl;
//          if(d==1)
//            std::cerr << "range<" << d << ">: " << k << " " << i << " " << j << std::endl;
//          if(d==2)
//            std::cerr << "range<" << d << ">: " << j << " " << k << " " << i << std::endl;

        using idxperm::pis;

        // second reconstructed position (j=2, between i=1 and i=2)
        const rng_t i_next = i + 2*dist;

        const int x[3] = {0, 1, 2};

        // helper references
        const arr_t u_0 = arr(pis<d>(i      - dist, j, k)),
                    u_1 = arr(pis<d>(i      + dist, j, k)),
                    u_2 = arr(pis<d>(i_next + dist, j, k));

        arr_t c_1 = c_j(pis<d>(i,      j, k)),
              c_2 = c_j(pis<d>(i_next, j, k)),
              d_1 = d_j(pis<d>(i,      j, k)),
              d_2 = d_j(pis<d>(i_next, j, k)),
              f_1 = f_j(pis<d>(i,      j, k)),
              f_2 = f_j(pis<d>(i_next, j, k));

        // Eqs. (5) and (6) Akinlabi et al.
        c_1 = (u_1 - u_0) / (2.) - d_1 * (u_2 - u_0) / (2.);
        c_2 = (u_2 - u_1) / (2.) - d_2 * (u_2 - u_0) / (2.);
        f_1 = (x[2] * u_0 - x[0] * u_1) / (2.) - d_1 * (x[2] * u_0 - x[0] * u_2) / (2.);
        f_2 = (x[2] * u_1 - x[0] * u_2) / (2.) - d_2 * (x[2] * u_0 - x[0] * u_2) / (2.);

        // new u, at j=1, between i=0 and i=1, result of w_j=1(u_i=1), Eq. (1) in Akinlabi et al.
        arr(pis<d>(i    , j, k)) = c_1 * 1 + d_1 * u_1 + f_1; 
        // new u, at j=2, between i=1 and i=2, result of w_j=2(u_i=1), Eq. (1) in Akinlabi et al.
        arr(pis<d>(i_next, j, k)) = c_2 * 1 + d_2 * u_1 + f_2;

        if(blitz::min(arr(pis<d>(i    , j, k))) < 0)
        {
          std::cerr << "negative reonstruction @ i: " << 
            " i:   " << i <<
            " j:   " << j <<
            " k:   " << k <<
            " arr: " << arr(pis<d>(i    , j, k)) <<
                  " u_0: " << u_0 <<
                  " u_1: " << u_1 <<
                  " u_2: " << u_2 <<
                  " c_1: " << c_1 <<
                  " d_1: " << d_1 <<
                  " f_1: " << f_1 <<
                  std::endl;
        }

        if(blitz::min(arr(pis<d>(i_next, j, k))) < 0)
        {
          std::cerr << "negative reonstruction @ i_next: " << 
            " i_next:" << i_next <<
            " j:   " << j <<
            " k:   " << k <<
            " arr: " << arr(pis<d>(i_next, j, k)) <<
                  " u_0: " << u_0 <<
                  " u_1: " << u_1 <<
                  " u_2: " << u_2 <<
                  " c_2: " << c_2 <<
                  " d_2: " << d_2 <<
                  " f_2: " << f_2 <<
                  std::endl;
        }

          // DEBUG: do interpolation, useful for testing if ranges are correct
//          arr(pis<d>(i, j, k)) = real_t(.5) * (u_0 + u_1);
//          arr(pis<d>(i_next, j, k)) = real_t(.5) * (u_1 + u_2);


// alternative solution following Scotti and Meneveau Physica D 1999
/*

          // helper references
          // NOTE: these are actually the resolved values only in the first iteration, in the subsequent iterations some of these are reconstructed values
          //       should we point to the resolved values in all iterations?? 
          const arr_t u_im1 = arr(pis<d>(i      - dist, j, k)), // resolved u_(i-1)
                      u_i   = arr(pis<d>(i      + dist, j, k)), // resolved u_(i)
                      u_ip1 = arr(pis<d>(i_next + dist, j, k)); // resolved u_(i+1)

          // c_j and f_j naming comes from the Akinlabi paper, TODO: rename to a_i b_i
          // NOTE: in multiple iterations of reconstruction, they only change due to random stretching parameters? 
          //       or the stretching parameter should be drawn from the distro once? then no need for c_j and f_j to have the same size as reconstructed fields...
          arr_t a_1 = this->c_j(pis<d>(i,      j, k)),
                a_2 = this->c_j(pis<d>(i_next, j, k)),
                b_1 = this->f_j(pis<d>(i,      j, k)),
                b_2 = this->f_j(pis<d>(i_next, j, k));

          // Eqs. (13) and (14) from Scotti Meneveau 1999 
          a_1 = u_i - u_im1 - d_1*(u_ip1 - u_im1);
          a_2 = u_ip1 - u_i - d_2*(u_ip1 - u_im1);
          b_1 = u_im1*(real_t(1) - d_1); 
          b_2 = u_i - d_2*u_im1; 

          // W(u_0(xi)) at xi=1/4, between u_i-1 and u_i
          // u_0(2*xi) = u_i-1 + 1/2 * (u_i+1 - u_i-1)
          // Eq. (6) and caption of fig. 1 Scotti 1999
          // NOTE: in multiple iterations, should xi be between the 3 resolved points? or should we consider the reconstructed points as new resolved points?
          arr(pis<d>(i,      j, k)) = d_1 * 
//                                      (u_im1 + 0.5 * (u_ip1 - u_im1)) + // u_0(2*xi)
                                      u_i + 
                                      0.5 * a_1 + b_1;                  // q_1(2*xi)
          // W(u_0(xi)) at xi=3/4, between u_i and u_i+1
          // u_0(2*xi-1) = u_i-1 + 1/2 * (u_i+1 - u_i-1)
          // NOTE: u_0(2*3/4-1) == u_0(2*1/4), so calculate it once in the previous step and store it...
          arr(pis<d>(i_next, j, k)) = d_2 * 
//                                      (u_im1 + 0.5 * (u_ip1 - u_im1)) + // u_0(2*xi-1)
                                      u_i + 
                                      0.5 * a_2 + b_2;                  // q_2(2*xi-1)
*/
          // Eq. (5) Scotti and Meneveau PRL 1997, a_1 and a_2 differ from Scotti Meneveau 1999 !
          // NOTE: in the 1999 paper, a_1 and a_2 are multipled by 2*xi and 2*xi-1, respectively
          //       in the 1997 paper, the are both multipled by xi
          /*
          a_1 = 2 * (u_im1 - u_i - d_1*(u_ip1 - u_im1));
          a_2 = 2 * (u_ip1 - u_i - d_2*(u_ip1 - u_im1));
          b_1 = u_im1*(real_t(1) - d_1); 
          b_2 = u_i - d_2*u_im1; 
          */
      }
    };
  };
};
