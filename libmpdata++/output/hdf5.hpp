/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-netCDF reader
 */

#pragma once

#include <libmpdata++/output/detail/output_timer.hpp>
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
    class hdf5 : public detail::output_timer<solver_t> // TODO: get rid of timer here!
    {
      using parent_t = detail::output_timer<solver_t>;

      //static_assert(parent_t::n_dims < 3, "only 1D and 2D output supported");

      const std::string outfile;
      std::unique_ptr<H5::H5File> hdfp;
      std::map<int, H5::DataSet> vars;
      std::map<int, H5::DataSet> dims;
      std::vector<hsize_t> t;

      // HDF types of host data
      const H5::FloatType 
        flttype_solver = 
	  sizeof(typename solver_t::real_t) == sizeof(long double) 
	    ? H5::PredType::NATIVE_LDOUBLE 
	    : sizeof(typename solver_t::real_t) == sizeof(double) 
	      ? H5::PredType::NATIVE_DOUBLE :
	      H5::PredType::NATIVE_FLOAT,
        flttype_output = H5::PredType::NATIVE_FLOAT;

      // TODO: consider saving to floats by default to conserve space?

      void start(const int nt)
      {
        try 
        {
          // turn off the default output printing
          //H5::Exception::dontPrint();

          // creating the file
          hdfp.reset(new H5::H5File(outfile, H5F_ACC_TRUNC));

          // creating the dimensions
          const int n_dims = parent_t::n_dims + 1; // +1 for time
          hsize_t shape[n_dims], limit[n_dims], chunk[n_dims];
          shape[0] = 0; 
          limit[0] = H5S_UNLIMITED;
          chunk[0] = 1;  // TODO: a better choice perhaps?
          for (int i = 1; i <= parent_t::n_dims; ++i) 
            shape[i] = limit[i] = chunk[i] = this->mem->span[i-1];

          // creating variables
	  // enabling chunking in order to use unlimited dimension
          {
            // some helpers
            const H5::StrType strtype(H5::PredType::C_S1, H5T_VARIABLE);
	    const hsize_t one = 1;
	    H5::DataSpace scalar(1, &one); // netcdf fails to read HDF5 scalars :(

            // creating the time variable
            {
              H5::DSetCreatPropList params;
              params.setChunk(1, chunk);

              H5::DataSpace space(1, shape, limit);
              dims[0] = (*hdfp).createDataSet("time", flttype_output, space, params);

herr_t status;
status = H5DSset_scale(dims[0].getId(), "time");
std::cerr << "status = " << status << std::endl;

              // it is meant to help Paraview guess which variable is time
              dims[0].createAttribute("unit", strtype, scalar).write(strtype, std::string("TODO since 0"));
            }

	    // creating the user-requested variables
            {
	      H5::DSetCreatPropList params;
	      params.setChunk(n_dims, chunk);
              params.setDeflate(5); // TODO: a better choice of algorithm or the parameter? // TODO: move such constant to the header

              H5::DataSpace space(n_dims, shape, limit);
	      for (const auto &v : outvars)
	      {
		vars[v.first] = (*hdfp).createDataSet(v.second.name, flttype_output, space, params);
		vars[v.first].createAttribute("unit", strtype, scalar).write(strtype, v.second.unit);
	      }
            }
          }
        }
        catch (const H5::Exception &e)
        {
          BOOST_THROW_EXCEPTION(
            error(e.getCDetailMsg()) 
              << boost::errinfo_api_function(e.getCFuncName())
              << boost::errinfo_file_name(outfile.c_str())
          );
        }
      }

      void stop()
      {
      }
 
      void record(const int var)
      {
        try
        {
	  switch (solver_t::n_dims)
	  {
	    case 2:
	    {
              // dimensions, offsets, etc
              t[var]++;
              const hsize_t 
                n0 = hsize_t(this->mem->span[0]),
                n1 = hsize_t(this->mem->span[1]),
                shape[3] = { t[var], n0, n1 };
              hsize_t
                count[3] = { 1,        1, n1 },
                start[3] = { t[var]-1, 0, 0  };

              vars[var].extend(shape);
              dims[0].extend(shape); // TODO: unnecesarily repeated for each var

              const H5::DataSpace column(3, count);
	      H5::DataSpace space = vars[var].getSpace();

	      // halos present -> data not contiguous -> looping over the major rank
	      for (int i = 0; i < this->mem->span[0]; ++i)
	      {
                start[1] = i;
                space.selectHyperslab( H5S_SELECT_SET, count, start);
		vars[var].write( &(this->mem->state(var)(i,0)), flttype_solver, column, space); 
	      }
	      break;
	    }
	    default: assert(false); // TODO: 1D and 3D versions
	  }
        }
        catch (const H5::Exception &e)
        {
          BOOST_THROW_EXCEPTION(
            error(e.getCDetailMsg()) 
              << boost::errinfo_api_function(e.getCFuncName())
              << boost::errinfo_file_name(outfile.c_str())
          );
        }
      }

      decltype(parent_t::params_t::outvars) outvars;

      public:

      struct params_t : parent_t::params_t 
      { 
	std::string outfile;
// TODO: pass adiitional info? (e.g. Thrust version for icicle)
      };

      // ctor
      hdf5(
	typename parent_t::ctor_args_t args,
	const params_t &p
      ) : parent_t(args, p), outfile(p.outfile), outvars(p.outvars), t(solver_t::n_eqs, 0)
      { }
    }; 
  }; // namespace output
}; // namespace libmpdataxx
