/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/stretching_parameter_data/Waclawczyk_LES/cdf_th_subsaturated.hpp>
#include <libmpdata++/formulae/stretching_parameter_data/Waclawczyk_LES/cdf_th_supersaturated.hpp>
#include <libmpdata++/formulae/stretching_parameter_data/Waclawczyk_LES/cdf_rv_subsaturated.hpp>
#include <libmpdata++/formulae/stretching_parameter_data/Waclawczyk_LES/cdf_rv_supersaturated.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace fractal
    {
      namespace stretch_params
      {
        // available distributions of the stretching parameter
        enum class d_distro_t {
          DNS_vel,
          LES_rv_subsaturated, 
          LES_rv_supersaturated, 
          LES_th_subsaturated, 
          LES_th_supersaturated
        };

        // stretching parameters for velocity from DNS
        // the function used is a fit to the E. Akinlabi data 
        // the data is in the stretching_parameter_data subdirectory
        namespace DNS_vel 
        {
          // CDF of stretching parameter d; assuming CDF = A d^B + C under condition that CDF(0.5)=0 and CDF(1)=1
          template<class real_t>
          static real_t CDF_of_d(const real_t &d)
          {
            const real_t B = -0.39091161; 
            return (pow(d,B) - pow(0.5,B)) / (1. - pow(0.5,B));
          }
          // inverse function: d of a CDF
          template<class real_t>
          static real_t d_of_CDF(const real_t &CDF)
          {
            const real_t B = -0.39091161; 
            return pow((1.-pow(0.5,B))*CDF + pow(0.5,B), 1./B);
          }
          template<class real_t>
          struct d_of_CDF_fctr_DNS
          {
            real_t operator()(const real_t &CDF) const
            {
              return CDF < 0 ? -d_of_CDF(-CDF) : d_of_CDF(CDF);
            }
            BZ_DECLARE_FUNCTOR(d_of_CDF_fctr_DNS);
          };
        };
        // stretching parameters for potential temperature in a LES of marine cumulus
        namespace LES_th_rv
        {
          template<class real_t>
          class d_of_CDF_fctr_LES
          {
            std::vector<std::pair<real_t, real_t>> d_CDF_vctr;

            public:

            d_of_CDF_fctr_LES(d_distro_t dd)
            {
              switch(dd)
              {
                case d_distro_t::LES_th_subsaturated:
                  d_CDF_th_subsaturated(d_CDF_vctr);
                  break;
                case d_distro_t::LES_th_supersaturated:
                  d_CDF_th_supersaturated(d_CDF_vctr);
                  break;
                case d_distro_t::LES_rv_subsaturated:
                  d_CDF_rv_subsaturated(d_CDF_vctr);
                  break;
                case d_distro_t::LES_rv_supersaturated:
                  d_CDF_rv_supersaturated(d_CDF_vctr);
                  break;
                default:
                  std::runtime_error("libmpdata++: invalid d_distro_t type in the constructor of d_of_CDF_fctr_LES");
              }
            };

            real_t operator()(const real_t &CDF) const
            {
              auto pos = std::lower_bound(d_CDF_vctr.begin(), d_CDF_vctr.end(), CDF, 
                [](const std::pair<real_t, real_t> &pair, real_t val)
                {
                  return pair.second < val;
                });
              return pos->first;
            }
            BZ_DECLARE_FUNCTOR(d_of_CDF_fctr_LES);
          };
        };
      };
    };
  };
};
