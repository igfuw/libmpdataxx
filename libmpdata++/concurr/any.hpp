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
    template <typename real_t, int n_dims>
    struct any
    {
      virtual 
      void advance(int) 
      { assert(false); throw; }  

      virtual 
      blitz::Array<real_t, n_dims> advectee(int eqn = 0)
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
      bool *panic_ptr() 
      { assert(false && "unimplemented!"); throw; }

      // dtor
      virtual ~any() {}
    };
  }
}
