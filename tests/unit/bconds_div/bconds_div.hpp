#pragma once
#include <libmpdata++/solvers/boussinesq.hpp>

template <class ct_params_t>
class bconds_div : public libmpdataxx::solvers::boussinesq<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::boussinesq<ct_params_t>;

  protected:
  
  const std::string error_str;

  void hook_post_step() final
  {
    auto gc_div = this->max_abs_vctr_div(this->mem->GC);

    if (gc_div > 2 * this->prs_tol)
    {
      if (this->rank == 0)
      {
        std::cout << "bconds: " << error_str
                  << " gc_div: " << gc_div << std::endl;
      }
      this->mem->barrier();
      throw std::runtime_error("");
    }
    
    parent_t::hook_post_step();
  }

  public:

  struct rt_params_t : parent_t::rt_params_t
  {
    std::string error_str;
  };

  bconds_div(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    error_str(p.error_str)
  {}
};
