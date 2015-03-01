/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-xdmf reader
 */


#pragma once

#include <libmpdata++/output/detail/output_common.hpp>
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
    class hdf5 : public detail::output_common<solver_t>
    {
      using parent_t = detail::output_common<solver_t>;

      protected:

      const std::string outdir;
      std::unique_ptr<H5::H5File> hdfp;
      std::map<int, H5::DataSet> vars;
      std::map<int, std::string> dim_names;

      // HDF types of host data
      const H5::FloatType
        flttype_solver =
	  sizeof(typename solver_t::real_t) == sizeof(long double)
	    ? H5::PredType::NATIVE_LDOUBLE
	    : sizeof(typename solver_t::real_t) == sizeof(double)
	      ? H5::PredType::NATIVE_DOUBLE :
	      H5::PredType::NATIVE_FLOAT,
        flttype_output = H5::PredType::NATIVE_FLOAT; // using floats not to waste disk space

      blitz::TinyVector<hsize_t, parent_t::n_dims> cshape, shape, chunk, count, offst;
      H5::DSetCreatPropList params;

      void start(const int nt)
      {
        {
          // creating the directory
          boost::filesystem::create_directory(outdir);

          // creating the coordinates file
          std::string dim_file = outdir + "/coord.h5";
std::cerr << "creating " << dim_file << "(rank=" << this->mem->rank() << ")" << std::endl;
          hdfp.reset(new H5::H5File(dim_file, H5F_ACC_TRUNC));
std::cerr << "... done (" << dim_file << ")." << std::endl;

          // creating the dimensions
          // x,y,z
          offst = 0;

          shape = this->mem->advectee().extent();
          
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
            // X, Y, Z
            for (int i = 0; i < parent_t::n_dims; ++i)
            {

              blitz::Array<typename solver_t::real_t, parent_t::n_dims> coord(cshape);
              std::string name;
              switch (i)
              {
                case 0 : coord = this->di * blitz::firstIndex();
                         name = "X";
                         dim_names[i] = name;
                         break;
                case 1 : coord = this->dj * blitz::secondIndex();
                         name = "Y";
                         dim_names[i] = name;
                         break;
                case 2 : coord = this->dk * blitz::thirdIndex();
                         name = "Z";
                         dim_names[i] = name;
                         break;
                default : break;
              }

              auto curr_dim = (*hdfp).createDataSet(name, flttype_output, H5::DataSpace(parent_t::n_dims, cshape.data()));

              H5::DataSpace dim_space = curr_dim.getSpace();
              dim_space.selectHyperslab(H5S_SELECT_SET, cshape.data(), offst.data());
              curr_dim.write(coord.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, cshape.data()), dim_space);
            }
            // T
            {
              const hsize_t 
                nt_out = nt / this->outfreq + 1, // incl. t=0
                zero = 0,
                one = 1;
              float dt = this->dt;

              blitz::Array<typename solver_t::real_t, 1> coord(nt_out);
              coord = this->outfreq * this->dt * blitz::firstIndex();

              auto curr_dim = (*hdfp).createDataSet("T", flttype_output, H5::DataSpace(1, &nt_out));

              H5::DataSpace dim_space = curr_dim.getSpace();
              dim_space.selectHyperslab(H5S_SELECT_SET, &nt_out, &zero);
              curr_dim.write(coord.data(), flttype_solver, H5::DataSpace(1, &nt_out), dim_space);

              curr_dim.createAttribute("dt", flttype_output, H5::DataSpace(1, &one)).write(flttype_output, &dt);
            }
          }
        }
      }

      std::string base_name()
      {
        std::stringstream ss;
        ss << "timestep" << std::setw(10) << std::setfill('0') << this->timestep;
        return ss.str();
      }
      
      std::string hdf_name()
      {
        // TODO: add option of .nc extension for Paraview sake ?
        return base_name() + ".h5";
      }

      void record_all()
      {
        // in concurrent setup only the first solver does output
        assert(this->mem->rank() == 0);
        //count[1] = 1; TODO

        // creating the timestep file
std::cerr << "creating "<< hdf_name() << " file ... (rank=" << this->mem->rank() << ")" << std::endl;
        hdfp.reset(new H5::H5File(outdir + "/" + hdf_name(), H5F_ACC_TRUNC));
std::cerr << "... done (" << hdf_name() << ")." << std::endl;

        {
	  for (const auto &v : this->outvars)
          {
            // creating the user-requested variables
            vars[v.first] = (*hdfp).createDataSet(
              v.second.name,
              flttype_output,
              H5::DataSpace(parent_t::n_dims, shape.data()),
              params
            );
	    // TODO: units attribute

            H5::DataSpace space = vars[v.first].getSpace();
            switch (int(solver_t::n_dims))
            {
              case 1:
              {
                space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                vars[v.first].write( &(this->mem->advectee(v.first)(0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), space);
                break;
              }
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
      }

      // data is assumed to be contiguous and in the same layout as hdf variable
      void record_aux(const std::string &name, typename solver_t::real_t *data)
      {
        assert(this->mem->rank() == 0);

        auto aux = (*hdfp).createDataSet(
          name,
          flttype_output,
          H5::DataSpace(parent_t::n_dims, shape.data()),
          params
        );

        auto space = aux.getSpace();
        offst = 0;
        space.selectHyperslab(H5S_SELECT_SET, shape.data(), offst.data());
        aux.write(data, flttype_solver, H5::DataSpace(parent_t::n_dims, shape.data()), space);
      }

      public:

      struct rt_params_t : parent_t::rt_params_t
      {
	std::string outdir;
// TODO: pass adiitional info? (command_line, library versions, ...) (-> output_common?)
      };

      // ctor
      hdf5(
	typename parent_t::ctor_args_t args,
	const rt_params_t &p
      ) : parent_t(args, p),
        outdir(p.outdir)
      { }
    };
  }; // namespace output
}; // namespace libmpdataxx
