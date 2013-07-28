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
      hsize_t t;
      int outfreq;

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

      void start(const int nt)
      {
        try 
        {
          // turn off the default output printing
          //H5::Exception::dontPrint(); // TODO: when done with boost exceptions set-up

          // creating the file
          hdfp.reset(new H5::H5File(outfile, H5F_ACC_TRUNC));

          // creating the dimensions
          hsize_t shape[hdf_dims], limit[hdf_dims], chunk[hdf_dims];
          // time
          shape[0] = 0; 
          limit[0] = H5S_UNLIMITED;
          chunk[0] = 1; 
          // x,y,z
          if (parent_t::n_dims == 2)
          {
            shape[1] = limit[1] = chunk[1] = this->mem->span[0];
            shape[2] = limit[2] = chunk[2] = 1;
	    shape[3] = limit[3] = chunk[3] = this->mem->span[1];
          }
          else assert(false && "TODO");

          // creating variables
	  // enabling chunking in order to use unlimited dimension
          {
            // some helpers
            const H5::StrType strtype(H5::PredType::C_S1);

            // creating the dimension variables
            {
              H5::DSetCreatPropList params;
              params.setChunk(1, chunk);

              dims[0] = (*hdfp).createDataSet("time", flttype_output, H5::DataSpace(1, shape, limit), params);

              // phony_dim_0 -> time
	      herr_t status;
	      status = H5DSset_scale(dims[0].getId(), "time");
              assert(status == 0);

              // it is meant to help Paraview guess which variable is time
              dims[0].createAttribute("axis", strtype, H5::DataSpace(H5S_SCALAR)).write(strtype, std::string("T"));
              dims[0].createAttribute("units", strtype, H5::DataSpace(H5S_SCALAR)).write(strtype, std::string("seconds"));

/* // TODO...
              // creating the D0 variable
              dims[1] = (*hdfp).createDataSet("D0", flttype_output, H5::DataSpace(1, shape+1));
	      status = H5DSset_scale(dims[1].getId(), "D0");
              assert(status == 0);
              for (hsize_t i = 0, one = 1; i < shape[1]; ++i)
              {
                float x = (i + .5) * this->dx;
		H5::DataSpace space = dims[1].getSpace();
		space.selectHyperslab(H5S_SELECT_SET, &one, &i);
                dims[1].write(&x, flttype_output, scalar, space);
              }
*/
            }

	    // creating the user-requested variables
            {
	      H5::DSetCreatPropList params;
	      params.setChunk(hdf_dims, chunk);
              params.setDeflate(5); // TODO: move such constant to the header

              H5::DataSpace space(hdf_dims, shape, limit);
	      for (const auto &v : outvars)
	      {
		vars[v.first] = (*hdfp).createDataSet(v.second.name, flttype_output, space, params);
		vars[v.first].createAttribute("units", strtype, H5::DataSpace(H5S_SCALAR)).write(strtype, v.second.unit);
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

      void record_all()
      {
	// dimensions, offsets, etc
	const hsize_t 
	  n0 = hsize_t(this->mem->span[0]),
	  n1 = hsize_t(this->mem->span[1]),
	  shape[hdf_dims] = { t+1, n0, 1, n1 };
	hsize_t
	  count[hdf_dims] = { 1,   1,  1, n1 },
	  start[hdf_dims] = { t,   0,  0, 0  };

        const H5::DataSpace column(hdf_dims, count);
        const H5::DataSpace scalar(1,      count);

        try
        {
	  dims[0].extend(shape);
	  {
	    float t_sec = t * outfreq * this->dt; // TODO: shouldn't it be a member of output common?
	    H5::DataSpace space = dims[0].getSpace();
	    space.selectHyperslab(H5S_SELECT_SET, count, start);
	    dims[0].write(&t_sec, flttype_output, scalar, space);
	  }

	  for (const auto &v : outvars) switch (solver_t::n_dims)
	  {
	    case 2:
	    {
              vars[v.first].extend(shape);

	      H5::DataSpace space = vars[v.first].getSpace();

	      // halos present -> data not contiguous -> looping over the major rank
	      for (int i = 0; i < this->mem->span[0]; ++i)
	      {
                start[1] = i;
                space.selectHyperslab( H5S_SELECT_SET, count, start);
		vars[v.first].write( &(this->mem->state(v.first)(i,0)), flttype_solver, column, space); 
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
        t++;
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
      ) : parent_t(args, p), 
        outfile(p.outfile), outvars(p.outvars), outfreq(p.outfreq), t(0) // TODO: all these should be members of output_common!
      { }
    }; 
  }; // namespace output
}; // namespace libmpdataxx
