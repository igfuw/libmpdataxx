#include <cmath>
#include <boost/math/constants/constants.hpp>

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

#include "test_def.hpp"

using namespace libmpdataxx;

int main()
{
  const bool var_dt = true;
  const T max_cfl = 0.8;

  std::vector<int> nys = {60, 120, 240, 480};

  if (FULL_SIM)
  {
    nys.push_back(960);
  }

  for (const auto ny : nys) 
  {
    {
      enum { opts = opts::nug | opts::iga | opts::fct};
      const int opts_iters = 2;
      test<var_dt, opts, opts_iters>("Mg2No", ny, max_cfl);
    }
    
    {
      enum { opts = opts::nug | opts::abs | opts::div_2nd | opts::div_3rd};
      const int opts_iters = 2;
      test<var_dt, opts, opts_iters>("Mp3", ny, max_cfl);
    }
    
    {
      enum { opts = opts::nug | opts::abs | opts::tot};
      const int opts_iters = 3;
      test<var_dt, opts, opts_iters>("Mp3cc", ny, max_cfl);
    }
    
    {
      enum { opts = opts::nug | opts::iga | opts::div_2nd | opts::div_3rd | opts::fct};
      const int opts_iters = 2;
      test<var_dt, opts, opts_iters>("Mg3No", ny, max_cfl);
    }
  }
}
