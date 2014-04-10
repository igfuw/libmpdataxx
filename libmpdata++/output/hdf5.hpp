/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-netCDF reader
 */

// TODO: rename hdf5_carthesian?

#pragma once

#include <libmpdata++/output/detail/output_common.hpp>
#include <libmpdata++/detail/error.hpp>

// the C++ HDF5 API
#include <H5Cpp.h>

// the C HDF5 API (for dimension scales)
#include "hdf5_hl.h"

namespace libmpdataxx
{
  namespace output
  {
    template <class solver_t>
    class hdf5 : public detail::output_common<solver_t> 
    {
      using parent_t = detail::output_common<solver_t>;

      //static_assert(parent_t::n_dims < 3, "only 1D and 2D output supported");

      const std::string outfile;
      std::unique_ptr<H5::H5File> hdfp;
      std::map<int, H5::DataSet> vars;
      std::map<int, H5::DataSet> dims;
      std::map<const std::string, H5::DataSet> aux;

      // HDF types of host data
      const H5::FloatType 
        flttype_solver = 
	  sizeof(typename solver_t::real_t) == sizeof(long double) 
	    ? H5::PredType::NATIVE_LDOUBLE 
	    : sizeof(typename solver_t::real_t) == sizeof(double) 
	      ? H5::PredType::NATIVE_DOUBLE :
	      H5::PredType::NATIVE_FLOAT,
        flttype_output = H5::PredType::NATIVE_FLOAT; // using floats not to waste disk space


      static const int hdf_dims = 4; // HDF dimensions (spatial + time)
      hsize_t shape[hdf_dims], limit[hdf_dims], chunk[hdf_dims], count[hdf_dims], offst[hdf_dims]; // TODO: std::arrays?
      H5::DSetCreatPropList params;

      void start(const int nt)
      {
        try 
        {
          // turn off the default output printing
          //H5::Exception::dontPrint(); // TODO: when done with boost exceptions set-up

          // creating the file
          hdfp.reset(new H5::H5File(outfile, H5F_ACC_TRUNC));

          // creating the dimensions
          // time
          shape[0] = 0; 
          limit[0] = H5S_UNLIMITED;
          chunk[0] = 1; 
          // x,y,z
          if (parent_t::n_dims == 2)
          {
            shape[1] = limit[1] = chunk[1] = this->mem->grid_size[0];
            shape[2] = limit[2] = chunk[2] = 1;
	    shape[3] = limit[3] = chunk[3] = this->mem->grid_size[1];

	    count[0] = count[1] = count[2] = 1;
            count[3] = shape[3];
	    offst[0] = offst[1] = offst[2] = offst[3] = 0;
          }
          else assert(false && "TODO");

	  params.setChunk(hdf_dims, chunk);
	  params.setDeflate(5); // TODO: move such constant to the header

          // creating variables
	  // enabling chunking in order to use unlimited dimension
          {
            // some helpers
            const H5::StrType strtype(H5::PredType::C_S1);

            // creating the dimension variables
            {
              H5::DSetCreatPropList time_params;
              time_params.setChunk(1, chunk);

              dims[0] = (*hdfp).createDataSet("time", flttype_output, H5::DataSpace(1, shape, limit), time_params);

              // phony_dim_0 -> time
	      herr_t status;
	      status = H5DSset_scale(dims[0].getId(), "time");
              assert(status == 0);

              // it is meant to help Paraview guess which variable is time
              dims[0].createAttribute("axis", strtype, H5::DataSpace(H5S_SCALAR)).write(strtype, std::string("T"));
              // TODO: units attribute

              // TODO...
              // creating the X,Y,Z variables
            }

	    // creating the user-requested variables
	    for (const auto &v : this->outvars)
	    {
	      vars[v.first] = (*hdfp).createDataSet(
                v.second.name, 
                flttype_output, 
                H5::DataSpace(hdf_dims, shape, limit), 
                params
              );
	      // TODO: units attribute
	    }
          }
        }
        catch (const H5::Exception &e) { handle(e); }
      }

      void record_all()
      {
        assert(this->mem->rank() == 0);

	// dimensions, offsets, etc
	shape[0] += 1;
        offst[0] = shape[0]-1;
        count[1] = 1;

        try
        {
	  dims[0].extend(shape);
	  {
	    float t_sec = this->timestep * this->dt; // TODO: shouldn't it be a member of output common?
	    H5::DataSpace time_space = dims[0].getSpace();
	    time_space.selectHyperslab(H5S_SELECT_SET, count, offst);
	    dims[0].write(&t_sec, flttype_output, H5::DataSpace(1, count), time_space);
	  }

	  for (const auto &v : this->outvars) switch (solver_t::n_dims)
	  {
	    case 2:
	    {
              vars[v.first].extend(shape);
	      H5::DataSpace space = vars[v.first].getSpace();

	      // halos present -> data not contiguous -> looping over the major rank
	      for (int i = 0; i < this->mem->grid_size[0]; ++i)
	      {
                offst[1] = i;
                space.selectHyperslab(H5S_SELECT_SET, count, offst);
		vars[v.first].write( &(this->mem->advectee(v.first)(i,0)), flttype_solver, H5::DataSpace(hdf_dims, count), space); 
	      }
	      break;
	    }
	    default: assert(false); // TODO: 1D and 3D versions
	  }
        }
        catch (const H5::Exception &e) { handle(e); }
      }

      void handle(const H5::Exception &e) 
      {
	BOOST_THROW_EXCEPTION(
	  error(e.getCDetailMsg()) 
	    << boost::errinfo_api_function(e.getCFuncName())
	    << boost::errinfo_file_name(outfile.c_str())
	);
      }

      protected:

      // auxiliary fields handling (e.g. diagnostic fields)
      void setup_aux(const std::string &name)
      {
        assert(this->mem->rank() == 0);
        try 
        {
	  aux[name] = (*hdfp).createDataSet(
	    name, 
	    flttype_output, 
	    H5::DataSpace(hdf_dims, shape, limit), 
	    params
	  );
        }
        catch (const H5::Exception &e) { handle(e); } 
      }

      // data is assumed to be contiguous and in the same layout as hdf variable
      void record_aux(const std::string &name, typename solver_t::real_t *data)
      {
        assert(this->mem->rank() == 0);
        try 
        {
          aux[name].extend(shape);
	  H5::DataSpace space = aux[name].getSpace();

	  offst[1] = 0;
          count[1] = shape[1];
	  space.selectHyperslab(H5S_SELECT_SET, count, offst);
	  aux[name].write(data, flttype_solver, H5::DataSpace(hdf_dims, count), space); 

        }
        catch (const H5::Exception &e) { handle(e); } 
      }

      public:

      struct rt_params_t : parent_t::rt_params_t 
      { 
	std::string outfile;
// TODO: pass adiitional info? (command_line, library versions, ...) (-> output_common?)
      };

      // ctor
      hdf5(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) : parent_t(args, p), 
        outfile(p.outfile)  
      { }
    }; 
  }; // namespace output
}; // namespace libmpdataxx
