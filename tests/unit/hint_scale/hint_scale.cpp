#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>

#include <cfenv> // for fesetround() - see TODO below 

using namespace libmpdataxx;

template <bool scale>
struct ct_params_t :ct_params_default_t
{
  using real_t = float;
  enum { n_dims = 2 };
  enum { n_eqns = 1 };
  enum { opts = opts::fct | opts::iga };

  static constexpr int hint_scale(const int &e) 
  {
    return scale ? 24 : 0;
  }
};

template <int bit>
bool test() 
{
  using solver_t = solvers::mpdata<ct_params_t<bit>>;
  typename solver_t::rt_params_t p;

  p.grid_size = {4,4};
  p.n_iters = 2;

  concurr::serial<solver_t, bcond::cyclic, bcond::cyclic, bcond::open, bcond::open> slv(p); 

  slv.advectee() =
              0, 0, 8.16638e+07, 8.16339e+07, NAN, NAN, 
    NAN, NAN, 0, 0, 7.8215e+07,  8.15976e+07, NAN, NAN,
    NAN, NAN, 0, 0, 8.16724e+07, 8.16339e+07, NAN, NAN, 
    NAN, NAN, 0, 0, 8.16638e+07, 8.16339e+07;
  auto max_init = max(slv.advectee());

  slv.advector(0) = 
              0.000572957855183631, 0.000286478927591815, -0.000286478927591815, -0.000572957855183631, NAN, NAN, 
    NAN, NAN, -0.00114591582678258, -0.000572957913391292, 0.000572957855183631, 0.00114591582678258,   NAN, NAN,
    NAN, NAN, 0.00057295779697597,  0.000286478898487985, -0.000286478869384154, -0.00057295779697597; 

  slv.advector(1) = 
              0,                     0,                     0,                     NAN, NAN,
    NAN, NAN, 0.000859436870086938,  0.00171887374017387,   0.000859436928294599,  NAN, NAN,
    NAN, NAN, -0.000859436753671616, -0.00171887350734323,  -0.000859436870086938, NAN, NAN,
    NAN, NAN, -1.70754757555791e-10, -3.41509515111582e-10, -1.70754757555791e-10; 

  slv.advance(1);
  std::cerr << slv.advectee();

  if (max(slv.advectee()) > max_init)
    throw std::runtime_error("values above initial max!");

  return min(slv.advectee()) >= 0;
}

int main()
{
  // TODO: find another sample data that does not require it
  fesetround(FE_DOWNWARD);
  
  // expecting fail
  if (test<false>()) 
    throw std::runtime_error("failed to detect expected values below zero!");

  // expecting success
  if (!test<true>())
    throw std::runtime_error("scale factors did not seem to have helped!");
}
