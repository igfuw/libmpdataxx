#include <advoocat/solvers/mpdata_2d.hpp>
#include <advoocat/solvers/solver_inhomo.hpp>

using namespace advoocat; // TODO: not here?

template <
  typename real_t, 
  int n_iters, 
  solvers::inhomo_e inhomo,
  int n_eqs
>
class cloud_common : public solvers::inhomo_solver<solvers::mpdata_2d<real_t, n_iters, n_eqs>, inhomo>
{
  using parent_t = solvers::inhomo_solver<solvers::mpdata_2d<real_t, n_iters, n_eqs>, inhomo>;

  protected:

  typename parent_t::arr_t rhod;

  public:

  struct params_t : parent_t::params_t 
  { 
    std::vector<real_t> rhod; // profile
  };

  // ctor
  cloud_common( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    rhod(args.mem->tmp[__FILE__][0][0])
  {
    assert(p.rhod.size() == this->j.last()+1);

    // initialising rhod array columnwise with data from the p.rhod profile
    for (int i = this->i.first(); i <= this->i.last(); ++i)
      for (int j = this->j.first(); j <= this->j.last(); ++j)
	rhod(i, j) = p.rhod[j];
  }  

  static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
  {
    using namespace advoocat::arakawa_c;
    parent_t::alloc(mem, nx, ny);
    const rng_t i(0, nx-1), j(0, ny-1);
    mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
    mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t(i^parent_t::halo, j^parent_t::halo)); // rhod
  }
};
