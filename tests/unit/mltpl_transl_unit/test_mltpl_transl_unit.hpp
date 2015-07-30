#pragma once

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

struct ct_params_t : libmpdataxx::ct_params_default_t
{
  enum { n_dims = 3 };
  enum { n_eqns = 3 };
  using real_t = double;
};

using slv_out_t = libmpdataxx::output::hdf5_xdmf<
  libmpdataxx::solvers::mpdata<
    ct_params_t
  >
>;
