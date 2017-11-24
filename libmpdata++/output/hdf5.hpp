/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 * @brief HDF5 output logic targetted at Paraview-xdmf reader
 */


#pragma once

#include <libmpdata++/output/detail/output_common.hpp>
#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp> // include the param2str maps and solver_family tags
#include <libmpdata++/solvers/boussinesq.hpp> // ditto

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

      using output_t = hdf5<solver_t>;

      std::unique_ptr<H5::H5File> hdfp;
      std::map<int, std::string> dim_names;
      const std::string const_name = "const.h5";
      std::string const_file;
      const hsize_t zero = 0, one = 1;

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

      void start(const typename parent_t::advance_arg_t nt)
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
          
          // save selected compile and runtime parameters, the choice depends on the solver family
//          record_params(*hdfp, typename parent_t::solver_family{});
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
                nt_out = nt / this->outfreq + 1; // incl. t=0
              float dt = this->dt;

              blitz::Array<typename solver_t::real_t, 1> coord(nt_out);
              coord = (this->var_dt ? this->outfreq : this->outfreq * this->dt) * blitz::firstIndex();

              auto curr_dim = (*hdfp).createDataSet("T", flttype_output, H5::DataSpace(1, &nt_out));

              H5::DataSpace dim_space = curr_dim.getSpace();
              dim_space.selectHyperslab(H5S_SELECT_SET, &nt_out, &zero);
              curr_dim.write(coord.data(), flttype_solver, H5::DataSpace(1, &nt_out), dim_space);
            }
            
            // G factor
            if (this->mem->G.get() != nullptr)
            {
              auto g_set = (*hdfp).createDataSet("G", flttype_output, sspace);
              record_dsc_helper(g_set, *this->mem->G);
            }
            
#if defined(USE_MPI)
            // close the const file with mpi access
            hdfp.release();
            // do record_params only on the first node
            if (this->mem->distmem.rank() == 0)
            {
              // open the const file with regular access
              hdfp.reset(new H5::H5File(const_file, H5F_ACC_TRUNC));
              // save selected compile and runtime parameters, the choice depends on the solver family
              record_params(*hdfp, typename parent_t::solver_family{});
            } 
            this->mem->distmem.barrier();
#else
            record_params(*hdfp, typename parent_t::solver_family{});
#endif
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

            record_dsc_helper(vars[v.first], this->out_data(v.first));
          }
        }
      }
      
      void record_dsc_helper(const H5::DataSet &dset, const typename solver_t::arr_t &arr)
      {
        H5::DataSpace space = dset.getSpace();

        switch (int(solver_t::n_dims))
        {
          case 1:
          {
            space.selectHyperslab(H5S_SELECT_SET, count.data(), offst.data());
            dset.write( &(arr(0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), space);
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
              dset.write( &(arr(i,0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), space);
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
                dset.write( &(arr(i,j,0)), flttype_solver, H5::DataSpace(parent_t::n_dims, count.data()), space);
              }
            }
            break;
          }
          default: assert(false);
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
      
      // for discontiguous array with halos
      void record_aux_dsc(const std::string &name, const typename solver_t::arr_t &arr)
      {
        assert(this->rank == 0);
        
        auto aux = (*hdfp).createDataSet(
          name,
          flttype_output,
          sspace,
          params
        );
        
        record_dsc_helper(aux, arr);
      }

      // has to be called after const file was created (i.e. after start())
      void record_aux_const(const std::string &name, typename solver_t::real_t data)
      {
        assert(this->rank == 0);
        float data_f(data);

        H5::H5File hdfcp(const_file, H5F_ACC_RDWR);; // reopen the const file
        hdfcp.openGroup("/").createAttribute(name, flttype_output, H5::DataSpace(1, &one)).write(flttype_output, &data_f);
      }

      // parameters saved for pure advection solvers
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_family_tag)
      {
        assert(this->rank == 0);
        hdfcp.createGroup("advection");
        const auto &group = hdfcp.openGroup("advection");
        {
          const auto opts_str = opts::opts_string(parent_t::ct_params_t_::opts);
          const auto type = H5::StrType(H5::PredType::C_S1, opts_str.size());
          group.createAttribute("opts", type, H5::DataSpace(1, &one)).write(type, opts_str.data());
        }
        {
          const auto type = flttype_solver;
          group.createAttribute("dt", type, H5::DataSpace(1, &one)).write(type, &this->dt);
          const auto names = std::vector<std::string>{"di", "dj", "dk"};
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            group.createAttribute(names[d], type, H5::DataSpace(1, &one)).write(type, &this->dijk[d]);
          }
        }
        {
          const auto type = H5::PredType::NATIVE_HBOOL;
          const auto data = parent_t::ct_params_t_::var_dt;
          group.createAttribute("var_dt", type, H5::DataSpace(1, &one)).write(type, &data);
        }
        if (parent_t::ct_params_t_::var_dt)
        {
          const auto type = flttype_solver;
          group.createAttribute("max_courant", type, H5::DataSpace(1, &one)).write(type, &this->max_courant);
        }
        {
          const auto type = H5::PredType::NATIVE_INT;
          const auto data = this->n_iters;
          group.createAttribute("n_iters", type, H5::DataSpace(1, &one)).write(type, &data);
        }
      }
      
      // as above but for solvers with rhs
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_rhs_family_tag)
      {
        record_params(hdfcp, typename solvers::mpdata_family_tag{});

        hdfcp.createGroup("rhs");
        const auto &group = hdfcp.openGroup("rhs");
        const auto scheme_str = solvers::scheme2string.at(static_cast<solvers::rhs_scheme_t>(parent_t::ct_params_t_::rhs_scheme));
        const auto type = H5::StrType(H5::PredType::C_S1, scheme_str.size());
        group.createAttribute("rhs_scheme", type, H5::DataSpace(1, &one)).write(type, scheme_str.data());
      }
      
      // as above but for solvers with velocities
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_rhs_vip_family_tag)
      {
        record_params(hdfcp, typename solvers::mpdata_rhs_family_tag{});
        
        hdfcp.createGroup("vip");
        const auto &group = hdfcp.openGroup("vip");
        const auto vab_str = solvers::vab2string.at(static_cast<solvers::vip_vab_t>(parent_t::ct_params_t_::vip_vab));
        const auto type = H5::StrType(H5::PredType::C_S1, vab_str.size());
        group.createAttribute("vip_abs", type, H5::DataSpace(1, &one)).write(type, vab_str.data());
      }
      
      // as above but for solvers with pressure equation
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_rhs_vip_prs_family_tag)
      {
        record_params(hdfcp, typename solvers::mpdata_rhs_vip_family_tag{});
        
        hdfcp.createGroup("prs");
        const auto &group = hdfcp.openGroup("prs");
        {
          const auto prs_scheme_str = solvers::prs2string.at(static_cast<solvers::prs_scheme_t>(parent_t::ct_params_t_::prs_scheme));
          const auto type = H5::StrType(H5::PredType::C_S1, prs_scheme_str.size());
          group.createAttribute("prs_scheme", type, H5::DataSpace(1, &one)).write(type, prs_scheme_str.data());
        }
        {
          const auto type = flttype_solver;
          group.createAttribute("prs_tol", type, H5::DataSpace(1, &one)).write(type, &this->prs_tol);
        }
      }
      
      // as above but for solvers with subgrid model (parameters common to all subgrid models)
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_rhs_vip_prs_sgs_family_tag)
      {
        record_params(hdfcp, typename solvers::mpdata_rhs_vip_prs_family_tag{});
        
        hdfcp.createGroup("sgs");
        const auto &group = hdfcp.openGroup("sgs");
        {
          const auto sgs_scheme_str = solvers::sgs2string.at(static_cast<solvers::sgs_scheme_t>(parent_t::ct_params_t_::sgs_scheme));
          const auto type = H5::StrType(H5::PredType::C_S1, sgs_scheme_str.size());
          group.createAttribute("sgs_scheme", type, H5::DataSpace(1, &one)).write(type, sgs_scheme_str.data());
        }
        {
          const auto sdiff_str = solvers::sdiff2string.at(static_cast<solvers::stress_diff_t>(parent_t::ct_params_t_::stress_diff));
          const auto type = H5::StrType(H5::PredType::C_S1, sdiff_str.size());
          group.createAttribute("stress_diff", type, H5::DataSpace(1, &one)).write(type, sdiff_str.data());
        }
        {
          const auto type = flttype_solver;
          group.createAttribute("cdrag", type, H5::DataSpace(1, &one)).write(type, &this->cdrag);
        }
      }
      
      // as above but for solvers with the dns subgrid model
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_rhs_vip_prs_sgs_dns_family_tag)
      {
        record_params(hdfcp, typename solvers::mpdata_rhs_vip_prs_sgs_family_tag{});
        
        const auto &group = hdfcp.openGroup("sgs");
        {
          const auto type = flttype_solver;
          group.createAttribute("eta", type, H5::DataSpace(1, &one)).write(type, &this->eta);
        }
      }
      
      // as above but for solvers with the smg subgrid model
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_rhs_vip_prs_sgs_smg_family_tag)
      {
        record_params(hdfcp, typename solvers::mpdata_rhs_vip_prs_sgs_family_tag{});
        
        const auto &group = hdfcp.openGroup("sgs");
        {
          const auto type = flttype_solver;
          group.createAttribute("smg_c", type, H5::DataSpace(1, &one)).write(type, &this->smg_c);
          group.createAttribute("c_m", type, H5::DataSpace(1, &one)).write(type, &this->c_m);
        }
      }
      
      // as above but for the boussinesq solver
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_boussinesq_family_tag)
      {
        record_params(hdfcp, typename std::conditional<static_cast<solvers::sgs_scheme_t>
                                                        (parent_t::ct_params_t_::sgs_scheme) == solvers::iles,
                                                       typename solvers::mpdata_rhs_vip_prs_family_tag,
                                                       typename solvers::mpdata_rhs_vip_prs_sgs_smg_family_tag>::type{});
        
        hdfcp.createGroup("boussinesq");
        const auto &group = hdfcp.openGroup("boussinesq");
        {
          const auto type = flttype_solver;
          group.createAttribute("g", type, H5::DataSpace(1, &one)).write(type, &this->g);
          group.createAttribute("Tht_ref", type, H5::DataSpace(1, &one)).write(type, &this->Tht_ref);
          group.createAttribute("hflux_const", type, H5::DataSpace(1, &one)).write(type, &this->hflux_const);
        }
        {
          auto dset = group.createDataSet("tht_e", flttype_output, H5::DataSpace(parent_t::n_dims, shape.data()));
          record_dsc_helper(dset, this->tht_e);
        }
      }

      // as above but for the boussinesq sgs solver
      void record_params(const H5::H5File &hdfcp, typename solvers::mpdata_boussinesq_sgs_family_tag)
      {
        record_params(hdfcp, typename solvers::mpdata_boussinesq_family_tag{});
        
        const auto &group = hdfcp.openGroup("boussinesq");
        {
          const auto type = flttype_solver;
          group.createAttribute("prandtl_num", type, H5::DataSpace(1, &one)).write(type, &this->prandtl_num);
        }
        {
          auto dset = group.createDataSet("mix_len", flttype_output, H5::DataSpace(parent_t::n_dims, shape.data()));
          record_dsc_helper(dset, this->mix_len);
        }
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
