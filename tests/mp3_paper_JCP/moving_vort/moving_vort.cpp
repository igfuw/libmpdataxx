#include <libmpdata++/solvers/mpdata.hpp>
#include "test_def.hpp"

using namespace libmpdataxx;

int main()
{
  std::vector<int> nys = {24, 48, 96, 192, 384, 768};

  const bool var_dt = true;
  const T max_cfl = 1.0;

  {
    enum { opts = opts::nug | opts::iga | opts::fct};
    const int opts_iters = 2;
    for (const auto ny : nys) test<var_dt, opts, opts_iters>("Mg2No", ny, max_cfl);
  }
  
  {
    enum { opts = opts::nug | opts::abs | opts::div_2nd | opts::div_3rd};
    const int opts_iters = 2;
    for (const auto ny : nys) test<var_dt, opts, opts_iters>("Mp3", ny, max_cfl);
  }
  
  {
    enum { opts = opts::nug | opts::abs | opts::dfl | opts::tot};
    const int opts_iters = 3;
    for (const auto ny : nys) test<var_dt, opts, opts_iters>("Mp3cc", ny, max_cfl);
  }
  
  {
    enum { opts = opts::nug | opts::iga | opts::div_2nd | opts::div_3rd | opts::fct};
    const int opts_iters = 2;
    for (const auto ny : nys) test<var_dt, opts, opts_iters>("Mg3No", ny, max_cfl);
  }
}
