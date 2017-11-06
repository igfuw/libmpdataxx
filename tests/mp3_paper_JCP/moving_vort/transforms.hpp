#pragma once

#include <boost/math/constants/constants.hpp>
#include <blitz/array.h>

using T = double;
constexpr T pi = boost::math::constants::pi<T>();

// forward transformations
struct xpf_t
{
  T x0, y0;
  T operator()(T x, T y) const
  { 
    T ret = atan2(cos(y) * sin(x - x0), cos(y) * sin(y0) * cos(x - x0) - cos(y0) * sin(y));
    return ret <= 0 ? ret + 2 * pi : ret;
  }
  BZ_DECLARE_FUNCTOR2(xpf_t);
};

struct ypf_t
{
  T x0, y0;
  T operator()(T x, T y) const
  { 
    return asin(sin(y) * sin(y0) + cos(y) * cos(y0) * cos(x - x0));
  }
  BZ_DECLARE_FUNCTOR2(ypf_t);
};

// inverse transformations
struct ixpf_t
{
  T x0, y0;
  T operator()(T x, T y) const
  { 

    T ret = x0 + atan2(cos(y) * sin(x), sin(y) * cos(y0) + cos(y) * cos(x) * sin(y0));
    return ret <= 0 ? ret + 2 * pi : ret;
  }
  BZ_DECLARE_FUNCTOR2(ixpf_t);
};

struct iypf_t
{
  T x0, y0;
  T operator()(T x, T y) const
  { 
    return asin(sin(y) * sin(y0) - cos(y) * cos(y0) * cos(x));
  }
  BZ_DECLARE_FUNCTOR2(iypf_t);
};
