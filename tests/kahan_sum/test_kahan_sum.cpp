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
  
  // that's why we might need a Kahan sum
  if (blitz::sum(a) - blitz::sum(a.reverse(0)) < 4e-8) throw std::exception();
  // that's to show it works
  if (blitz::kahan_sum(a) - blitz::kahan_sum(a.reverse(0)) > 4e-8) throw std::exception();
}
