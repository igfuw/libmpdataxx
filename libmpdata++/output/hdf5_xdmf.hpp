/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-xdmf reader
 */


#pragma once

#include <libmpdata++/output/hdf5.hpp>
#include <libmpdata++/output/detail/xdmf_writer.hpp>

#include <boost/filesystem.hpp>

// the C++ HDF5 API
#include <H5Cpp.h>

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
namespace libmpdataxx
{
  namespace output
  {
    template <class solver_t>
    class hdf5_xdmf : public hdf5<solver_t>
    {
      using parent_t = hdf5<solver_t>;

      static_assert(parent_t::n_dims > 1, "only 2D and 3D output supported");
      
      std::vector<std::string> timesteps;
      std::vector<std::string> attr_names;
      //xdmf writer
      detail::xdmf_writer<parent_t::n_dims> xdmfw;


      void start(const typename parent_t::advance_arg_t nt)
      {
        parent_t::start(nt);

        // get variable names for xdmf writer setup
        for (const auto &v : this->outvars)
        {
          attr_names.push_back(v.second.name);
        }
        
        if (this->mem->G.get() != nullptr) xdmfw.add_const_attribute("G", this->const_name, this->cshape);

        xdmfw.setup(this->const_name, this->dim_names, attr_names, this->cshape);
      }

      void record_all()
      {
        // write xdmf markup
        std::string xmf_name = this->base_name() + ".xmf";
        xdmfw.write(this->outdir + "/" + xmf_name, this->hdf_name(), this->record_time);

        // save the xmf filename for temporal write
        timesteps.push_back(xmf_name);
        // write temporal xmf
        xdmfw.write_temporal(this->outdir + "/temp.xmf", timesteps);

        parent_t::record_all();
      }

      public:

      struct rt_params_t : parent_t::rt_params_t
      {
        std::vector<std::string> outvars_aux;
      };

      // ctor
      hdf5_xdmf(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) : parent_t(args, p), attr_names(p.outvars_aux)
      {}
    };
  } // namespace output
} // namespace libmpdataxx
