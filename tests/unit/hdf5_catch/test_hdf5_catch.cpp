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
    enum { n_dims = 2 };
    enum { n_eqns = 1 };
    enum { rhs_scheme = solvers::euler_a };
  };

  // solver & output choice
  using solver_t = output::hdf5<solvers::mpdata_rhs<ct_params_t>>;

  // run-time parameters
  solver_t::rt_params_t p; 

  p.grid_size = { 10, 10 };
  p.dt = 1;
  p.outfile = "-/-"; // <- intentionaly invalid

  // instantiation
  concurr::serial<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p);

  // disabling automatic printing of exceptions
  auto error_stack = H5Eget_current_stack();
  H5Eset_auto(error_stack, NULL, NULL);

  // integration
  try
  {
    run.advance(1); 
  } 
  catch (H5::Exception)
  {
    std::cerr << "error caught!" << std::endl;
    exit(EXIT_SUCCESS);
  }
  exit(EXIT_FAILURE);
};
