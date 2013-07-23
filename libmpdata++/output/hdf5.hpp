/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/output/detail/output_timer.hpp>
#include <libmpdata++/detail/error.hpp>

#include <H5Cpp.h>

namespace libmpdataxx
{
  namespace output
  {
    template <class solver_t>
    class hdf5 : public detail::output_timer<solver_t>
    {
      using parent_t = detail::output_timer<solver_t>;

      //static_assert(parent_t::n_dims < 3, "only 1D and 2D output supported");

      //std::unique_ptr<Gnuplot> gp;
      H5::H5File file;

      void start(const int nt)
      {
        //gp.reset(new Gnuplot());

      }

      void stop()
      {
        //gp.reset();
      }
 
      void record(const int var)
      {
      }

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
      ) : parent_t(args, p)
      {
        try 
        {
          // turn off the default output printing
          H5::Exception::dontPrint();
          file.openFile(p.outfile, H5F_ACC_TRUNC);
          //file.createDataSet();
        }
        catch (const H5::Exception &e)
        {
          BOOST_THROW_EXCEPTION(
            error(e.getCDetailMsg()) 
              << boost::errinfo_api_function(e.getCFuncName())
              << boost::errinfo_file_name(p.outfile.c_str())
          );
        }
      }
    }; 
  }; // namespace output
}; // namespace libmpdataxx
