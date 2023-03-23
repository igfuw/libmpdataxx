/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-xdmf reader
 */


#pragma once

#include <libmpdata++/output/hdf5.hpp>

namespace libmpdataxx
{
  namespace output
  {
    // specialization for solvers with grid refinemenet
    template <class solver_t>
    class hdf5<solver_t,
      typename std::enable_if_t<(int)solver_t::ct_params_t::fractal_recon != (int)0>
    > : public hdf5_common<solver_t>
    {
      protected:

      using parent_t = hdf5_common<solver_t>;
      using parent_t::parent_t;
      using output_t = hdf5<solver_t>;

      blitz::TinyVector<hsize_t, parent_t::n_dims> shape_ref, cshape_ref, chunk_ref, offst_ref;
      H5::DataSpace sspace_ref, cspace_ref;

      std::map<int, std::string> dim_names_ref;

      void start(const typename parent_t::advance_arg_t nt) override
      {
        parent_t::start(nt);

        offst_ref = 0;
        for (int d = 0; d < parent_t::n_dims; ++d)
          shape_ref[d] = this->mem->distmem.grid_size_ref[d];

        chunk_ref  = shape_ref;
        cshape_ref = shape_ref + 1;
        sspace_ref = H5::DataSpace(parent_t::n_dims, shape_ref.data());
        cspace_ref = H5::DataSpace(parent_t::n_dims, cshape_ref.data());

#if defined(USE_MPI)
        if (this->mem->distmem.size() > 1)
        {
          shape_ref[0] = this->mem->distmem.rank() < this->mem->distmem.size() - 1 ? 
            this->mem->grid_size_ref[0].length()-1 : // -1 to make ranges nonoverlapping
            this->mem->grid_size_ref[0].length();
          cshape_ref[0] = shape_ref[0];

          if (this->mem->distmem.rank() == this->mem->distmem.size() - 1)
            cshape_ref[0] += 1;

          offst_ref[0] = this->mem->grid_size_ref[0].first();
          chunk_ref[0] = ( (typename solver_t::real_t) (this->mem->distmem.grid_size_ref[0])) / this->mem->distmem.size() + 0.5 ;
        }
#endif

        // refined X, Y, Z, TODO: very similar to X, Y, Z
        for (int i = 0; i < parent_t::n_dims; ++i)
        {
          blitz::Array<typename solver_t::real_t, parent_t::n_dims> coord(this->cshape_ref);
#if defined(USE_MPI)
          coord.reindexSelf(this->offst_ref);
#endif
          std::string name;
          switch (i)
          {
            case 0 : coord = this->di / 2. + this->di / this->mem->n_ref * (blitz::firstIndex() - .5);
                     name = "X refined";
                     dim_names_ref[i] = name;
                     break;
            case 1 : coord = this->dj / 2. + this->dj / this->mem->n_ref * (blitz::secondIndex() - .5);
                     name = "Y refined";
                     dim_names_ref[i] = name;
                     break;
            case 2 : coord = this->dk / 2. + this->dk / this->mem->n_ref * (blitz::thirdIndex() - .5);
                     name = "Z refined";
                     dim_names_ref[i] = name;
                     break;
            default : break;
          }
          auto curr_dim = (*(this->hdfp)).createDataSet(name, this->flttype_output, this->cspace_ref);

          H5::DataSpace dim_space = curr_dim.getSpace();
          dim_space.selectHyperslab(H5S_SELECT_SET, this->cshape_ref.data(), this->offst_ref.data());
          curr_dim.write(coord.data(), this->flttype_solver, H5::DataSpace(parent_t::n_dims, this->cshape_ref.data()), dim_space, this->dxpl_id);
        }
      }

      void record_dsc_ref_helper(const H5::DataSet &dset, const typename solver_t::arr_t &arr)
      {
        parent_t::record_dsc_helper(dset, arr, this->mem->grid_size_ref, shape_ref, offst_ref);
      }

      void record_aux_ref_hlpr(const std::string &name, typename solver_t::real_t *data, H5::H5File &hdf)
      {
        this->params.setChunk(parent_t::n_dims, chunk_ref.data());

        record_aux_hlpr(name, data, hdf5, sspace_ref, shape_ref, offst_ref);

        // revert to default chunk
        this->params.setChunk(parent_t::n_dims, this->chunk.data());
      }

      void record_aux_dsc_ref_hlpr(const std::string &name, const typename solver_t::arr_t &arr, H5::H5File hdf)
      {
        parent_t::record_aux_dsc_hlpr(name, arr, hdf, false, 
          record_dsc_ref_helper, 
          chunk_ref, 
          sspace_ref
        );
      }




      void record_aux_dsc_refined(const std::string &name, const typename solver_t::arr_t &arr)
      {
        record_aux_dsc_ref_hlpr(name, arr, *hdfp);
      }

      void record_aux_refined(const std::string &name, typename solver_t::real_t *data)
      {
        record_aux_ref_hlpr(name, data, *hdfp);
      }

    };

  } // namespace output
} // namespace libmpdataxx
