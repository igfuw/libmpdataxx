/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "concurrent_1d/test_concurrent_1d.cpp"
 */

#include <libmpdata++/concurr/openmp.hpp>
#include <libmpdata++/concurr/boost_thread.hpp>
#include <libmpdata++/concurr/cxx11_thread.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/concurr/threads.hpp>

#include <libmpdata++/solvers/mpdata.hpp>

#include <libmpdata++/output/hdf5.hpp>

int main()
{
#if defined(USE_MPI)
  // we will instantiate many solvers, so we have to init mpi manually, 
  // because solvers will not know should they finalize mpi upon destruction
  MPI::Init_thread(MPI_THREAD_MULTIPLE);
#endif
  using namespace libmpdataxx;

  std::cerr << "OpenMP: ";
#if defined(_OPENMP)
  std::cerr << "on" << std::endl;
#else 
  std::cerr << "off" << std::endl;
#endif

  struct ct_params_t : ct_params_default_t
  { 
    using real_t = double; 
    enum { n_dims = 1 }; 
    enum { n_eqns = 1 };
  };

  const int  
    nx = 100, // if it is less than the number of cores, the default setting will fail
    nt = 1000;
   
  // OpenMP
  std::cerr << "OpenMP run" << std::endl;
  {
    using slv_t = solvers::mpdata<ct_params_t>;
    using slv_out_t = output::hdf5<slv_t>;
    using run_t = concurr::openmp<
      slv_out_t, 
      bcond::cyclic, bcond::cyclic
    >;
    typename slv_out_t::rt_params_t p;
    p.grid_size = { nx };
    p.outdir = "output";
    run_t run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }

  // C++11 threads
  std::cerr << "C++11 threads run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::cxx11_thread<solver_t, bcond::cyclic, bcond::cyclic> run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }

  // Boost.Thread
  std::cerr << "Boost.Thread run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::boost_thread<solver_t, bcond::cyclic, bcond::cyclic> run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }

  // trheads (i.e. auto)
  std::cerr << "threads run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }

  // serial
  std::cerr << "serial run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::serial<solver_t, bcond::cyclic, bcond::cyclic> run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }
#if defined(USE_MPI)
  MPI::Finalize();
#endif
}
