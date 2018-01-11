#include <libmpdata++/solvers/mpdata.hpp>
using namespace libmpdataxx;

#include "../common/convergence.hpp"
#include "../../mp3_paper_2018_JCP/moving_vort/test_def.hpp"

int main()
{
#if defined(USE_MPI)
  // we will instantiate many solvers, so we have to init mpi manually, 
  // because solvers will not know should they finalize mpi upon destruction
  MPI::Init_thread(MPI_THREAD_MULTIPLE);
#endif
  const bool var_dt = true;
  const T max_cfl = 0.99;
  
  {
    enum { opts = opts::nug | opts::tot};
    const int opts_iters = 3;
    convergence(test<var_dt, opts, opts_iters>, "nug_tot_i3", max_cfl);
  }
  
  {
    enum { opts = opts::nug | opts::iga | opts::tot | opts::fct};
    const int opts_iters = 2;
    convergence(test<var_dt, opts, opts_iters>, "nug_iga_tot_fct_i2", max_cfl);
  }
#if defined(USE_MPI)
  MPI::Finalize();
#endif

}

