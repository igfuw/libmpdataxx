#include "advoocat/blitz.hpp"

int main()
{
  blitz::Array<float, 1> a(10);
  a = 1e-8;
  a(0) = 1e8;
  std::cerr  
    << a
    << std::setprecision(20)
    << "blitz::sum = " << blitz::sum(a) << std::endl
    << " kahan_sum = " << blitz::kahan_sum(a) << std::endl;
}
