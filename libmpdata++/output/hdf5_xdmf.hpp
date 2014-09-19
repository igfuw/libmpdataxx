/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-xdmf reader
 */


#pragma once

#include <libmpdata++/output/detail/output_common.hpp>
#include <libmpdata++/detail/error.hpp>
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
    class hdf5_xdmf : public detail::output_common<solver_t>
    {
      using parent_t = detail::output_common<solver_t>;

      //static_assert(parent_t::n_dims < 3, "only 1D and 2D output supported");

      const std::string outdir;
      std::unique_ptr<H5::H5File> hdfp;
      std::map<int, H5::DataSet> vars;
      std::map<int, H5::DataSet> dims;
      std::map<const std::string, H5::DataSet> aux;
      std::map<int, std::string> dim_names;
      std::vector<std::string> timesteps;

      // HDF types of host data
      const H5::FloatType
        flttype_solver =
	  sizeof(typename solver_t::real_t) == sizeof(long double)
	    ? H5::PredType::NATIVE_LDOUBLE
	    : sizeof(typename solver_t::real_t) == sizeof(double)
	      ? H5::PredType::NATIVE_DOUBLE :
	      H5::PredType::NATIVE_FLOAT,
        flttype_output = H5::PredType::NATIVE_FLOAT; // using floats not to waste disk space

      blitz::TinyVector<hsize_t, parent_t::n_dims> cshape, shape, limit, chunk, count, offst;
      H5::DSetCreatPropList params;

      //xdmf writer
      detail::xdmf_writer<parent_t::n_dims> xdmfw;

      void start(const int nt)
      {
        try
        {
          // turn off the default output printing
          //H5::Exception::dontPrint(); // TODO: when done with boost exceptions set-up

          // creating the directory
          boost::filesystem::create_directory(outdir);

          // creating the coordinates file
          std::string dim_file = outdir + "/coord.h5";
          hdfp.reset(new H5::H5File(dim_file, H5F_ACC_TRUNC));

          // creating the dimensions
          // x,y,z
          offst = 0;

          limit = shape = this->mem->advectee().extent();
          
          chunk = 1;
          // change chunk size along the last dimension
          *(chunk.end() - 1) = *(shape.end() - 1);

          count = 1;
          // see above
          *(count.end() - 1) = *(shape.end() - 1);

          // there is one more coordinate than cell index in each dimension
          cshape = shape + 1;

	  params.setChunk(parent_t::n_dims, chunk.data());
	  params.setDeflate(5); // TODO: move such constant to the header

          // creating variables
          {
            for (int i = 0; i < parent_t::n_dims; ++i)
            {

              blitz::Array<typename solver_t::real_t, parent_t::n_dims> coord(cshape);
              std::string name;
              switch (i)
              {
                case 0 : coord = blitz::firstIndex();
                         name = "X";
                         dim_names[i] = name;
                         break;
                case 1 : coord = blitz::secondIndex();
                         name = "Y";
                         dim_names[i] = name;
                         break;
                case 2 : coord = blitz::thirdIndex();
                         name = "Z";
                         dim_names[i] = name;
                         break;
                default : break;
              }

              dims[i] = (*hdfp).createDataSet(name, flttype_output, H5::DataSpace(parent_t::n_dims, cshape.data()));

              H5::DataSpace dim_space = dims[i].getSpace();
              dim_space.selectHyperslab(H5S_SELECT_SET, cshape.data(), offst.data());
              dims[i].write(coord.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, cshape.data()), dim_space);
            }
          }

          // get variable names for xdmf writer setup
          std::vector<std::string> attr_names;
          for (const auto &v : this->outvars)
          {
            attr_names.push_back(v.second.name);
          }

          xdmfw.setup("coord.h5", dim_names, attr_names, cshape);
        }
        catch (const H5::Exception &e) { handle(e); }
      }

      void record_all()
      {
        // in concurrent setup only the first solver does output
        assert(this->mem->rank() == 0);
        //count[1] = 1; TODO

        std::stringstream base_name;
        base_name << "timestep" << std::setw(10) << std::setfill('0') << this->timestep;

        std::string hdf_name = base_name.str() + ".h5";
        // creating the timestep file
        hdfp.reset(new H5::H5File(outdir + "/" + hdf_name, H5F_ACC_TRUNC));

        // write xdmf markup
        std::string xmf_name = base_name.str()+ ".xmf";
        xdmfw.write(outdir + "/" + xmf_name, hdf_name, this->timestep);

        // save the xmf filename for temporal write
        timesteps.push_back(xmf_name);
        // write temporal xmf
        xdmfw.write_temporal(outdir + "/temp.xmf", timesteps);

        try
        {
	  for (const auto &v : this->outvars)
          {
            // creating the user-requested variables
            vars[v.first] = (*hdfp).createDataSet(
              v.second.name,
              flttype_output,
              H5::DataSpace(parent_t::n_dims, shape.data(), limit.data()),
              params
            );
	    // TODO: units attribute

            H5::DataSpace space = vars[v.first].getSpace();
            switch (int(solver_t::n_dims))
            {
              case 2:
              {
                // halos present -> data not contiguous -> looping over the major rank
                for (int i = 0; i < this->mem->grid_size[0]; ++i)
                {
                  offst[0] = i;
                  space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                  vars[v.first].write( &(this->mem->advectee(v.first)(i,0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), space);
                }
                break;
              }
              case 3:
              {
                // halos present -> data not contiguous -> looping over the major rank
                for (int i = 0; i < this->mem->grid_size[0]; ++i)
                {
                  for (int j = 0; j < this->mem->grid_size[1]; ++j)
                  {
                    offst[0] = i;
                    offst[1] = j;
                    space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                    vars[v.first].write( &(this->mem->advectee(v.first)(i,j,0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), space);
                  }
                }
                break;
              }
              default: assert(false); // TODO: 1D version
            }
          }
        }
        catch (const H5::Exception &e) { handle(e); }
      }

      void handle(const H5::Exception &e)
      {
	BOOST_THROW_EXCEPTION(
	  error(e.getCDetailMsg())
	    << boost::errinfo_api_function(e.getCFuncName())
	    << boost::errinfo_file_name(outdir.c_str())
	);
      }

      protected:
      // TODO
      // auxiliary fields handling (e.g. diagnostic fields)
      //void setup_aux(const std::string &name)
      //{
      //  assert(this->mem->rank() == 0);
      //  try
      //  {
      //    aux[name] = (*hdfp).createDataSet(
      //      name,
      //      flttype_output,
      //      H5::DataSpace(hdf_dims, shape, limit),
      //      params
      //    );
      //  }
      //  catch (const H5::Exception &e) { handle(e); }
      //}

      //// data is assumed to be contiguous and in the same layout as hdf variable
      //void record_aux(const std::string &name, typename solver_t::real_t *data)
      //{
      //  assert(this->mem->rank() == 0);
      //  try
      //  {
      //    aux[name].extend(shape);
      //    H5::DataSpace space = aux[name].getSpace();

      //    offst[1] = 0;
      //    count[1] = shape[1];
      //    space.selectHyperslab(H5S_SELECT_SET, count, offst);
      //    aux[name].write(data, flttype_solver, H5::DataSpace(hdf_dims, count), space);

      //  }
      //  catch (const H5::Exception &e) { handle(e); }
      //}

      public:

      struct rt_params_t : parent_t::rt_params_t
      {
	std::string outdir;
// TODO: pass adiitional info? (command_line, library versions, ...) (-> output_common?)
      };

      // ctor
      hdf5_xdmf(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) : parent_t(args, p),
        outdir(p.outdir)
      { }
    };
  }; // namespace output
}; // namespace libmpdataxx
