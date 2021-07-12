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

      blitz::TinyVector<hsize_t, parent_t::n_dims> cshape, shape, chunk, srfcshape, srfcchunk, offst;
      H5::DSetCreatPropList params;

      H5::DataSpace sspace, cspace, srfcspace;
#if defined(USE_MPI)
      hid_t fapl_id;
#endif
      hid_t dxpl_id;

      void start(const typename parent_t::advance_arg_t nt)
      {
        const_file = this->outdir + "/" + const_name;

        if (this->mem->distmem.rank() == 0)
        {
          // creating the directory
          boost::filesystem::create_directory(this->outdir);
        }


        {
          // creating the const file
#if defined(USE_MPI)
          this->mem->distmem.barrier();
          H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
          H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
#endif
          hdfp.reset(new H5::H5File(const_file, H5F_ACC_TRUNC
#if defined(USE_MPI)
            , H5P_DEFAULT, fapl_id
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

          chunk = shape;

          // there is one more coordinate than cell index in each dimension
          cshape = shape + 1;

          srfcshape = shape;
          *(srfcshape.end()-1) = 1;
          sspace = H5::DataSpace(parent_t::n_dims, shape.data());
          srfcspace = H5::DataSpace(parent_t::n_dims, srfcshape.data());
          cspace = H5::DataSpace(parent_t::n_dims, cshape.data());

#if defined(USE_MPI)
          if (this->mem->distmem.size() > 1)
          {
            shape[0] = this->mem->grid_size[0].length();
            cshape[0] = this->mem->grid_size[0].length();

            if (this->mem->distmem.rank() == this->mem->distmem.size() - 1)
              cshape[0] += 1;

            offst[0] = this->mem->grid_size[0].first();

            // chunk size has to be common to all processes !
            // TODO: something better ?
            chunk[0] = ( (typename solver_t::real_t) (this->mem->distmem.grid_size[0])) / this->mem->distmem.size() + 0.5 ;
          }
#endif

          srfcshape = shape;
          *(srfcshape.end()-1) = 1;

          srfcchunk = chunk;
          *(srfcchunk.end()-1) = 1;

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
              curr_dim.write(coord.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, cshape.data()), dim_space, dxpl_id);
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
              curr_dim.write(coord.data(), flttype_solver, H5::DataSpace(1, &nt_out), dim_space, dxpl_id);
            }

            // G factor
            if (this->mem->G.get() != nullptr)
            {
              auto g_set = (*hdfp).createDataSet("G", flttype_output, sspace);
              record_dsc_helper(g_set, *this->mem->G);
            }

            // save selected compile and runtime parameters, the choice depends on the solver family
            record_params(*hdfp, typename parent_t::solver_family{});
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
            , H5P_DEFAULT, fapl_id
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

      void record_dsc_srfc_helper(const H5::DataSet &dset, const typename solver_t::arr_t &arr)
      {
        H5::DataSpace space = dset.getSpace();
        space.selectHyperslab(H5S_SELECT_SET, srfcshape.data(), offst.data());
        // TODO: some permutation of grid_size instead of the switch
        blitz::Range zro(0,0);

        switch (int(solver_t::n_dims))
        {
          case 1:
          {
            typename solver_t::arr_t contiguous_arr(1);
            contiguous_arr = arr(0); // create a copy that is contiguous
            dset.write(contiguous_arr.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, srfcshape.data()), space, dxpl_id);
            break;
          }
          case 2:
          {
            typename solver_t::arr_t contiguous_arr(this->mem->grid_size[0], zro);
            contiguous_arr = arr(this->mem->grid_size[0], zro); // create a copy that is contiguous
            dset.write(contiguous_arr.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, srfcshape.data()), space, dxpl_id);
            break;
          }
          case 3:
          {
            // create a copy that is contiguous and has the C-style (kji) storage order as required by HDF5
            typename solver_t::arr_t contiguous_arr(this->mem->grid_size[0], this->mem->grid_size[1], zro); 
            contiguous_arr = arr(this->mem->grid_size[0], this->mem->grid_size[1], zro);

            dset.write(contiguous_arr.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, srfcshape.data()), space, dxpl_id);
            break;
          }
          default: assert(false);
        };
      }

      void record_dsc_helper(const H5::DataSet &dset, const typename solver_t::arr_t &arr)
      {
        H5::DataSpace space = dset.getSpace();
        space.selectHyperslab(H5S_SELECT_SET, shape.data(), offst.data());
        // TODO: some permutation of grid_size instead of the switch

        switch (int(solver_t::n_dims))
        {
          case 1:
          {
            typename solver_t::arr_t contiguous_arr = arr(this->mem->grid_size[0]).copy(); // create a copy that is contiguous
            dset.write(contiguous_arr.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, shape.data()), space, dxpl_id);
            break;
          }
          case 2:
          {
            typename solver_t::arr_t contiguous_arr = arr(this->mem->grid_size[0], this->mem->grid_size[1]).copy(); // create a copy that is contiguous
            dset.write(contiguous_arr.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, shape.data()), space, dxpl_id);
            break;
          }
          case 3:
          {
            // create a copy that is contiguous and has the C-style (kji) storage order as required by HDF5
            typename solver_t::arr_t contiguous_arr(this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]); 
            contiguous_arr = arr(this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]);

            dset.write(contiguous_arr.data(), flttype_solver, H5::DataSpace(parent_t::n_dims, shape.data()), space, dxpl_id);
            break;
          }
          default: assert(false);
        };
      }

      // data is assumed to be contiguous and in the same layout as hdf variable and in the C-style storage order
      void record_aux_hlpr(const std::string &name, typename solver_t::real_t *data, H5::H5File hdf)
      {
        assert(this->rank == 0);

        auto aux = hdf.createDataSet(
          name,
          flttype_output,
          sspace,
          params
        );

        auto space = aux.getSpace();
        space.selectHyperslab(H5S_SELECT_SET, shape.data(), offst.data());
        aux.write(data, flttype_solver, H5::DataSpace(parent_t::n_dims, shape.data()), space, dxpl_id);
      }

      // data is assumed to be contiguous and in the same layout as hdf variable and in the C-style storage order
      void record_aux(const std::string &name, typename solver_t::real_t *data)
      {
        record_aux_hlpr(name, data, *hdfp);
      }

      // for discontiguous array with halos
      void record_aux_dsc_hlpr(const std::string &name, const typename solver_t::arr_t &arr, H5::H5File hdf, bool srfc = false)
      {
        assert(this->rank == 0);

        if(srfc)
          params.setChunk(parent_t::n_dims, srfcchunk.data());

        auto aux = hdf.createDataSet(
          name,
          flttype_output,
          srfc ? srfcspace : sspace,
          params
        );

        if(srfc)
        {
          // revert to default chunk
          params.setChunk(parent_t::n_dims, chunk.data());
          record_dsc_srfc_helper(aux, arr);
        }
        else
          record_dsc_helper(aux, arr);
      }

      void record_aux_dsc(const std::string &name, const typename solver_t::arr_t &arr, bool srfc = false)
      {
        record_aux_dsc_hlpr(name, arr, *hdfp, srfc);
      }


      void record_scalar_hlpr(const std::string &name, const std::string &group_name, typename solver_t::real_t data, H5::H5File hdf)
      {
        assert(this->rank == 0);
        float data_f(data);

        H5::Group group;
        // open a group if it exists, create it if it doesn't exist
        // based on: https://stackoverflow.com/questions/35668056/test-group-existence-in-hdf5-c
        // note: pre Hdf5-1.10, H5Lexists returns 0 for root group, hence we check directly if it is the root group
        // (https://support.hdfgroup.org/HDF5/doc/RM/RM_H5L.html#Link-Exists)
        if (group_name == "/" || H5Lexists(hdf.getId(), group_name.c_str(), H5P_DEFAULT) > 0)
          group = hdf.openGroup(group_name);
        else
          group = hdf.createGroup(group_name);

        group.createAttribute(name, flttype_output, H5::DataSpace(1, &one)).write(flttype_output, &data_f);
      }

      // has to be called after const file was created (i.e. after start())
      void record_aux_const(const std::string &name, typename solver_t::real_t *data)
      {
        H5::H5File hdfcp(const_file, H5F_ACC_RDWR); // reopen the const file
        record_aux_hlpr(name, data, hdfcp);
      }

      // has to be called after const file was created (i.e. after start())
      void record_aux_dsc_const(const std::string &name,  const typename solver_t::arr_t &arr)
      {
        H5::H5File hdfcp(const_file, H5F_ACC_RDWR); // reopen the const file
        record_aux_dsc_hlpr(name, arr, hdfcp);
      }

      void record_aux_const(const std::string &name, const std::string &group_name, typename solver_t::real_t data)
      {
        H5::H5File hdfcp(const_file, H5F_ACC_RDWR
#if defined(USE_MPI)
          , H5P_DEFAULT, fapl_id
#endif
        ); // reopen the const file
        record_scalar_hlpr(name, group_name, data, hdfcp);
      }

      void record_aux_const(const std::string &name, typename solver_t::real_t data)
      {
        record_aux_const(name, "/", data);
      }

      void record_aux_scalar(const std::string &name, const std::string &group_name, typename solver_t::real_t data)
      {
        record_scalar_hlpr(name, group_name, data, *hdfp);
      }

      void record_aux_scalar(const std::string &name, typename solver_t::real_t data)
      {
        record_aux_scalar(name, "/", data);
      }

      // see above, also assumes that z is the last dimension
      void record_prof_const(const std::string &name, typename solver_t::real_t *data)
      {
        assert(this->rank == 0);

        H5::H5File hdfcp(const_file, H5F_ACC_RDWR
#if defined(USE_MPI)
          , H5P_DEFAULT, fapl_id
#endif
        ); // reopen the const file

        auto aux = hdfcp.createDataSet(
          name,
          flttype_output,
          H5::DataSpace(1, &shape[parent_t::n_dims - 1])
        );

#if defined(USE_MPI)
        if (this->mem->distmem.rank() == 0)
#endif
        {
          auto space = aux.getSpace();
          space.selectHyperslab(H5S_SELECT_SET, &shape[parent_t::n_dims - 1], &offst[parent_t::n_dims - 1]);
          aux.write(data, flttype_solver, H5::DataSpace(1, &shape[parent_t::n_dims - 1]), space);
        }
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
          auto dset = group.createDataSet("tht_e", flttype_output, sspace);
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
          auto dset = group.createDataSet("mix_len", flttype_output, sspace);
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
        fapl_id = H5Pcreate(H5P_FILE_ACCESS);
#endif
        dxpl_id = H5Pcreate(H5P_DATASET_XFER);

        // TODO: clean it up - it should not be here
        // overrding the default from output_common
        if (this->outvars.size() == 1 && parent_t::n_eqns == 1)
          this->outvars[0].name = "psi";
      }

      // dtor
      virtual ~hdf5()
      {
        H5Pclose(dxpl_id);
#if defined(USE_MPI)
        H5Pclose(fapl_id);
#endif
      }
    };
  } // namespace output
} // namespace libmpdataxx
