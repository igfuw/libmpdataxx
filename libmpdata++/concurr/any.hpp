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
//<listing-1>
      void advance(real_t) 
//</listing-1>
      { assert(false); throw; }  

      virtual 
//<listing-2>
      blitz::Array<real_t, n_dims> advectee(int eqn = 0)
//</listing-2>
      { assert(false); throw; }

      virtual 
//<listing-3>
      blitz::Array<real_t, n_dims> advector(int dim = 0) 
//</listing-3>
      { assert(false); throw; }

      virtual 
//<listing-4>
      blitz::Array<real_t, n_dims> g_factor() 
//</listing-4>
      { assert(false); throw; }
     
      virtual 
      blitz::Array<real_t, n_dims> vab_coefficient() 
      { assert(false); throw; }
      
      virtual 
      blitz::Array<real_t, n_dims> vab_relaxed_state(int d = 0) 
      { assert(false); throw; }

      virtual 
//<listing-5>
      bool *panic_ptr() 
//</listing-5>
      { assert(false && "unimplemented!"); throw; }

      // dtor
      virtual ~any() {}
    };
  };
};
