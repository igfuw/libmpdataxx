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
      protected:

      using output_t = hdf5_xdmf<solver_t>;
      using parent_t = hdf5<solver_t>;

      static_assert(parent_t::n_dims > 1, "only 2D and 3D output supported");

      std::vector<std::string> timesteps;
      //xdmf writer
      detail::xdmf_writer<parent_t::n_dims> xdmfw;

      void start(const typename parent_t::advance_arg_t nt)
      {
        parent_t::start(nt);

#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
        {
          // get variable names for xdmf writer setup
          std::vector<std::string> attr_names;
          for (const auto &v : this->outvars)
          {
            attr_names.push_back(v.second.name);
          }

          if (this->mem->G.get() != nullptr) xdmfw.add_const_attribute("G", this->const_name, this->mem->distmem.grid_size.data());

          xdmfw.setup(this->const_name, this->dim_names, attr_names, this->mem->distmem.grid_size.data());
        }
      }

      void write_xmfs()
      {
#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
        {
          // write xdmf markup
          std::string xmf_name = this->base_name() + ".xmf";
          xdmfw.write(this->outdir + "/" + xmf_name, this->hdf_name(), this->record_time);

          // save the xmf filename for temporal write
          timesteps.push_back(xmf_name);
          // write temporal xmf
          xdmfw.write_temporal(this->outdir + "/temp.xmf", timesteps);
        }
      }

      void record_all()
      {
        write_xmfs();
        parent_t::record_all();
      }

      void record_aux(const std::string &name, typename solver_t::real_t *data)
      {
#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
          xdmfw.add_attribute(name, this->hdf_name(), this->mem->distmem.grid_size.data());
        parent_t::record_aux(name, data);
      }

      void record_aux_dsc(const std::string &name, const typename solver_t::arr_t &arr, bool srfc = false)
      {
        auto shape = this->mem->distmem.grid_size;
        if(srfc) shape.at(parent_t::n_dims-1) = 1;
#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
          xdmfw.add_attribute(name, this->hdf_name(), shape.data());
        parent_t::record_aux_dsc(name, arr, srfc);
      }

      public:

      // ctor
      hdf5_xdmf(
        typename parent_t::ctor_args_t args,
        const typename parent_t::rt_params_t &p
      ) : parent_t(args, p)
      {}
    };
  } // namespace output
} // namespace libmpdataxx
