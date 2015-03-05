// unit test for HDF5 exception passing through libmpdata++
//
// author[s]: Sylwester Arabas
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include <libmpdata++/solvers/mpdata_rhs.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/hdf5.hpp>

using namespace libmpdataxx; 

int main() 
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqns = 1 };
    enum { rhs_scheme = solvers::euler_a };
  };

  // solver & output choice
  using solver_t = output::hdf5<solvers::mpdata<ct_params_t>>;

  // run-time parameters
  solver_t::rt_params_t p; 

  p.grid_size = { 10 };
  p.outdir = boost::filesystem::unique_path().native();

  // instantiation
  concurr::serial<
    solver_t, 
    bcond::cyclic, bcond::cyclic
  > run(p);

  // integration
  try
  {
    run.advector() = 0;
    run.advectee() = 0;
    run.advance(1); 
    boost::filesystem::rename(boost::filesystem::path(p.outdir), boost::filesystem::path(p.outdir + ".off"));
    run.advance(1);  // this will fail and trow a HDF exception
  } 
  catch (H5::Exception)
  {
    std::cerr << "error caught!" << std::endl;
    exit(EXIT_SUCCESS);
  }
  exit(EXIT_FAILURE);
};
