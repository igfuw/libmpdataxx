/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-xdmf reader
 */


#pragma once

#include <libmpdata++/output/detail/output_common.hpp>

#include <boost/filesystem.hpp>

// the C++ HDF5 API
#include <H5Cpp.h>

#if defined(USE_MPI) && !defined(H5_HAVE_PARALLEL)
#  error "MPI enabled in libmpdata++ but not in HDF5"
#endif

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

      std::unique_ptr<H5::H5File> hdfp;
      std::map<int, std::string> dim_names;
      const std::string const_name = "const.h5";

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

      H5::DataSpace sspace, cspace;
#if defined(USE_MPI)
      hid_t plist_id;
#endif

      void start(const int nt)
      {
        std::string const_file = this->outdir + "/" + const_name;

        if (this->mem->distmem.rank() == 0)
        {
          // creating the directory
          boost::filesystem::create_directory(this->outdir);
        }

 
        {
          // creating the const file
#if defined(USE_MPI)
          this->mem->distmem.barrier();
          H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
          // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // TODO: check!
#endif
          hdfp.reset(new H5::H5File(const_file, H5F_ACC_TRUNC
#if defined(USE_MPI)
            , H5P_DEFAULT, plist_id
#endif
          ));
        }

        {

          // creating the dimensions
          // x,y,z
          offst = 0;

          for (int d = 0; d < parent_t::n_dims; ++d)
	    shape[d] = this->mem->distmem.grid_size[d];

          chunk = 1;
          // change chunk size along the last dimension
          *(chunk.end() - 1) = *(shape.end() - 1);

          count = 1;
          // see above
          *(count.end() - 1) = *(shape.end() - 1);

          // there is one more coordinate than cell index in each dimension
          cshape = shape + 1;

          sspace = H5::DataSpace(parent_t::n_dims, shape.data());
          cspace = H5::DataSpace(parent_t::n_dims, cshape.data());

#if defined(USE_MPI)
          if (this->mem->distmem.size() > 1)
          {
	    shape[0] = this->mem->grid_size[0].length();
	    cshape[0] = this->mem->grid_size[0].length();

	    if (this->mem->distmem.rank() == this->mem->distmem.size() - 1) 
              cshape[0] += 1;

            offst[0] = this->mem->grid_size[0].first();

            if (parent_t::n_dims == 1)
            {
              // chunk size has to be common to all processes !
              // TODO: something better ?
              chunk[0] = this->mem->distmem.grid_size[0] / this->mem->distmem.size();
              count[0] = shape[0];
            }
          }
#endif

	  params.setChunk(parent_t::n_dims, chunk.data());
#if !defined(USE_MPI)
	  params.setDeflate(5); // TODO: move such constant to the header
#endif

          // creating variables
          {
            // X, Y, Z
            for (int i = 0; i < parent_t::n_dims; ++i)
            {

              blitz::Array<typename solver_t::real_t, parent_t::n_dims> coord(cshape);
#if defined(USE_MPI)
              coord.reindexSelf(offst);
#endif
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

              auto curr_dim = (*hdfp).createDataSet(name, flttype_output, cspace);

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

            // G factor
            if (this->mem->G.get() != nullptr)
            {
              auto g_set = (*hdfp).createDataSet("G", flttype_output, sspace);
              H5::DataSpace g_space = g_set.getSpace();
              switch (int(solver_t::n_dims))
              {
                case 1:
                {
                  g_space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                  g_set.write( &((*this->mem->G)(0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), g_space);
                  break;
                }
                case 2:
                {
                  // halos present -> data not contiguous -> looping over the major rank
		  assert(this->mem->grid_size[0].stride() == 1);
		  for (auto i  = this->mem->grid_size[0].first(); 
			    i <= this->mem->grid_size[0].last(); 
			    i += this->mem->grid_size[0].stride()
                  ) {
                    offst[0] = i;
                    g_space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                    g_set.write( &((*this->mem->G)(i,0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), g_space);
                  }
                  break;
                }
                case 3:
                {
                  // halos present -> data not contiguous -> looping over the major rank
		  assert(this->mem->grid_size[0].stride() == 1);
		  for (auto i  = this->mem->grid_size[0].first(); 
			    i <= this->mem->grid_size[0].last(); 
			    i += this->mem->grid_size[0].stride()
                  ) {
		    assert(this->mem->grid_size[1].stride() == 1);
		    for (auto j  = this->mem->grid_size[1].first(); 
			      j <= this->mem->grid_size[1].last(); 
			      j += this->mem->grid_size[1].stride()
                    ) {
                      offst[0] = i;
                      offst[1] = j;
                      g_space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                      g_set.write( &((*this->mem->G)(i,j,0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), g_space);
                    }
                  }
                  break;
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
        assert(this->rank == 0);
        //count[1] = 1; TODO

	// creating the timestep file
	hdfp.reset(new H5::H5File(this->outdir + "/" + hdf_name(), H5F_ACC_TRUNC
#if defined(USE_MPI)
            , H5P_DEFAULT, plist_id
#endif
	));

        {
          std::map<int, H5::DataSet> vars;

	  for (const auto &v : this->outvars)
	  {
	    // creating the user-requested variables
	    vars[v.first] = (*hdfp).createDataSet(
	      v.second.name,
	      flttype_output,
	      sspace,
	      params
	    );
	    // TODO: units attribute

            H5::DataSpace space = vars[v.first].getSpace();
            switch (int(solver_t::n_dims))
            {
              case 1:
              {
                space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                vars[v.first].write( 
                  (this->mem->advectee(v.first).dataFirst()), 
                  flttype_solver, 
                  H5::DataSpace(parent_t::n_dims, count.data()), 
                  space
                );
                break;
              }
              case 2:
              {
                // halos present -> data not contiguous -> looping over the major rank
	        assert(this->mem->grid_size[0].stride() == 1);
	        for (auto i  = this->mem->grid_size[0].first(); 
                          i <= this->mem->grid_size[0].last(); 
                          i += this->mem->grid_size[0].stride()
                ) {
                  offst[0] = i;
                  space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
                  vars[v.first].write( &(this->mem->advectee(v.first)(i,0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), space);
                }
                break;
              }
              case 3:
              {
                // halos present -> data not contiguous -> looping over the major rank
	        assert(this->mem->grid_size[0].stride() == 1);
	        for (auto i  = this->mem->grid_size[0].first(); 
                          i <= this->mem->grid_size[0].last(); 
                          i += this->mem->grid_size[0].stride()
                ) {
	          assert(this->mem->grid_size[1].stride() == 1);
	          for (auto j  = this->mem->grid_size[1].first(); 
                            j <= this->mem->grid_size[1].last(); 
                            j += this->mem->grid_size[1].stride()
                  ) {
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
        assert(this->rank == 0);

        auto aux = (*hdfp).createDataSet(
          name,
          flttype_output,
          sspace,
          params
        );

        auto space = aux.getSpace();
        offst = 0;
        space.selectHyperslab(H5S_SELECT_SET, shape.data(), offst.data());
        aux.write(data, flttype_solver, H5::DataSpace(parent_t::n_dims, shape.data()), space);
      }

      public:

      // ctor
      hdf5(
	typename parent_t::ctor_args_t args,
	const typename parent_t::rt_params_t &p
      ) : parent_t(args, p)
      {
#if defined(USE_MPI)
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
#endif
        //plist_id = H5Pcreate(H5Pcreate(H5P_DATASET_XFER); // check!

        // TODO: clean it up - it should not be here
        // overrding the default from output_common
        if (this->outvars.size() == 1 && parent_t::n_eqns == 1)
          this->outvars[0].name = "psi";
      }

      // dtor
      virtual ~hdf5()
      {
#if defined(USE_MPI)
        H5Pclose(plist_id);
#endif
      }
    };
  } // namespace output
} // namespace libmpdataxx
