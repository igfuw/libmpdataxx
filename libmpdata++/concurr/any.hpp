// unit test for distributed memory facility
//
// author[s]: Sylwester Arabas
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/blitz.hpp>

namespace libmpdataxx
{
  namespace concurr
  {
    template <typename real_t, int n_dims, typename advance_arg_t = int>
    struct any
    {
      virtual 
      void advance(advance_arg_t) 
      { assert(false); throw; }  

      virtual 
      blitz::Array<real_t, n_dims> advectee(int eqn = 0)
      { assert(false); throw; }

      virtual 
      const blitz::Array<real_t, n_dims> advectee_global(int eqn = 0)
      { assert(false); throw; }

      virtual 
      void advectee_global_set(const blitz::Array<real_t, n_dims>, int eqn = 0)
      { assert(false); throw; }

      virtual 
      blitz::Array<real_t, n_dims> advector(int dim = 0) 
      { assert(false); throw; }

      virtual 
      blitz::Array<real_t, n_dims> g_factor() 
      { assert(false); throw; }
     
      virtual 
      blitz::Array<real_t, n_dims> vab_coefficient() 
      { assert(false); throw; }
      
      virtual 
      blitz::Array<real_t, n_dims> vab_relaxed_state(int d = 0) 
      { assert(false); throw; }
      
      virtual 
      blitz::Array<real_t, n_dims> sclr_array(const std::string &name, int n = 0)
      { assert(false); throw; }

      virtual 
      bool *panic_ptr() 
      { assert(false && "unimplemented!"); throw; }
      
      virtual 
      const real_t time() const
      { assert(false); throw; }
      
      // minimum of an advectee, mpi-aware
      virtual 
      const real_t min(int eqn = 0) const
      { assert(false); throw; }
      
      // maximum of an advectee, mpi-aware
      virtual 
      const real_t max(int eqn = 0) const
      { assert(false); throw; }

      // dtor
      virtual ~any() {}
    };
  }
}
