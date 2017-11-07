#include <libmpdata++/blitz.hpp>
#include <stdexcept>

template <typename real_t>
void main_tmpl()
{
  ////////////////////////////////////////////////////////////////
  // Kahan sum reduction test
  if (sizeof(real_t) == sizeof(float)) 
  {
    blitz::Array<real_t, 1> a(10), rev(10);
    a = 1e-8;
    a(0) = 1e8;
    std::cerr << std::setprecision(20);
    std::cerr  
      << a
      << "      blitz::sum = " << blitz::sum(a) << std::endl
      << "       kahan_sum = " << blitz::kahan_sum(a) << std::endl;
    rev = a.reverse(0);
    std::cerr  
      << rev
      << "rev + blitz::sum = " << blitz::sum(rev) << std::endl
      << "rev +  kahan_sum = " << blitz::kahan_sum(rev) << std::endl;
    
    // expected without Kahan sum
    auto expected_diff = 1.4e-8;

    // that's why we might need a Kahan sum
    if (std::abs(blitz::sum(a) - blitz::sum(rev)) < expected_diff) 
      throw std::runtime_error("A");
    // that's to show it works
    if (std::abs(blitz::kahan_sum(a) - blitz::kahan_sum(rev)) > expected_diff)
      throw std::runtime_error("B");
  }
}

int main()
{
  main_tmpl<float>();
  main_tmpl<double>();
  main_tmpl<long double>();
}
