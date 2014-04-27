/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/blitz.hpp>

// local definition of a Kahan's sum reduction
// (http://en.wikipedia.org/wiki/Kahan_summation_algorithm)
BZ_NAMESPACE(blitz)
  template<typename P_sourcetype, typename P_resulttype = BZ_SUMTYPE(P_sourcetype)>
  class ReduceKahanSum 
  {
    public:

    typedef P_sourcetype T_sourcetype;
    typedef P_resulttype T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceKahanSum() { } 

#pragma GCC push_options
#pragma GCC optimize ("O3") // assuming -Ofast could optimise out the algorithm
    bool operator()(const T_sourcetype& x, const int=0) const 
    { 
      T_resulttype y, t;
      y = x - c_;
      t = sum_ + y;
      c_ = (t - sum_) - y;
      sum_ = t;
      return true;
    }   
#pragma GCC pop_options

    T_resulttype result(const int) const { return sum_; }

    void reset() const 
    { 
      sum_ = c_ = zero(T_resulttype()); 
    }
 
    static const char* name() { return "sum"; }
 
    protected:

    mutable T_resulttype sum_, c_;
  };
  BZ_DECL_ARRAY_PARTIAL_REDUCE(kahan_sum, ReduceKahanSum)
  BZ_DECL_ARRAY_FULL_REDUCE(kahan_sum, ReduceKahanSum)
BZ_NAMESPACE_END
