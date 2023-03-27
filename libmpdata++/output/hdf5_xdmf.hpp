/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-xdmf reader
 */


#pragma once

#include <libmpdata++/output/hdf5_ref.hpp>
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
    class hdf5_xdmf_common : public hdf5<solver_t>
    {
      protected:

      using parent_t = hdf5<solver_t>;

      static_assert(parent_t::n_dims > 1, "only 2D and 3D output supported");

      std::vector<std::string> timesteps;
      //xdmf writer, additional one for refined data (separate .xmf files describing the refined grid, TODO: use two grids in one file? Paraview AMR dataset could help?)
      detail::xdmf_writer<parent_t::n_dims> xdmfw;

      void start(const typename parent_t::advance_arg_t nt) override
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

          xdmfw.setup(    this->const_name, this->dim_names,     attr_names, this->mem->distmem.grid_size.data()    );
        }
      }

      virtual void write_xmfs()
      {
#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
        {
          // write xdmf markup
          std::string xmf_name = this->base_name() + ".xmf";
          xdmfw.write(this->outdir + "/" + xmf_name, this->hdf_name(), this->record_time);

          // save the xmf filename for temporal write
          timesteps.push_back(this->base_name());
          // write temporal xmf
          xdmfw.write_temporal(this->outdir + "/temp.xmf", timesteps, ".xmf");
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
      hdf5_xdmf_common(
        typename parent_t::ctor_args_t args,
        const typename parent_t::rt_params_t &p
      ) : parent_t(args, p)
      {}
    };

    //template<class solver_t, class enableif = void>
    //class hdf5_xdmf;
           
    // for solvers without grid refinemenet
    //template <class solver_t>
    //class hdf5_xdmf<solver_t,
    //  typename std::enable_if_t<~libmpdataxx::solvers::detail::slvr_with_frac_recn_v<solver_t::ct_params_t>>
    //>
    
    template<class solver_t, class enableif = void>
    class hdf5_xdmf
    : public hdf5_xdmf_common<solver_t>
    {
      protected:
      using output_t = hdf5_xdmf<solver_t>;
      using parent_t = hdf5_xdmf_common<solver_t>;

      public:
      using parent_t::parent_t;
    };

    // specialization for solvers with grid refinemenet
    template <class solver_t>
    class hdf5_xdmf<solver_t,
      typename std::enable_if_t<(int)solver_t::ct_params_t_::fractal_recon != (int)0>
    > : public hdf5_xdmf_common<solver_t>
    {
      protected:
      using output_t = hdf5_xdmf<solver_t>;
      using parent_t = hdf5_xdmf_common<solver_t>;

      detail::xdmf_writer<parent_t::n_dims> xdmfw_ref;

      void start(const typename parent_t::advance_arg_t nt) override
      {
        parent_t::start(nt);

#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
        {
          xdmfw_ref.setup(this->const_name, this->dim_names_ref, {},         this->mem->distmem.grid_size_ref.data());
        }
      }

      void write_xmfs() override
      {
        parent_t::write_xmfs();
#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
        {
          std::string xmf_ref_name = this->base_name() + "_ref.xmf";
          xdmfw_ref.write(this->outdir + "/" + xmf_ref_name, this->hdf_name(), this->record_time);

          xdmfw_ref.write_temporal(this->outdir + "/temp_ref.xmf", this->timesteps, "_ref.xmf"); 
        }
      }

      void record_aux_refined(const std::string &name, typename solver_t::real_t *data)
      {
#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
          xdmfw_ref.add_attribute(name, this->hdf_name(), this->mem->distmem.grid_size_ref.data());
        parent_t::record_aux_refined(name, data);
      }

      void record_aux_dsc_refined(const std::string &name, const typename solver_t::arr_t &arr)
      {
        auto shape = this->mem->distmem.grid_size_ref;
#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
          xdmfw_ref.add_attribute(name, this->hdf_name(), shape.data());
        parent_t::record_aux_dsc_refined(name, arr);
      }

      public:
      using parent_t::parent_t;
    };


  } // namespace output
} // namespace libmpdataxx
