#include <libmpdata++/blitz.hpp>

template <typename real_t>
void main_tmpl()
{
  ////////////////////////////////////////////////////////////////
  // Kahan sum reduction test
  if (sizeof(real_t) == sizeof(float)) 
  {
    blitz::Array<real_t, 1> a(10);
    a = 1e-8;
    a(0) = 1e8;
    std::cerr  
      << a
      << std::setprecision(20)
      << "blitz::sum = " << blitz::sum(a) << std::endl
      << " kahan_sum = " << blitz::kahan_sum(a) << std::endl;
    
    // that's why we might need a Kahan sum
    if (blitz::sum(a) - blitz::sum(a.reverse(0)) < 4e-8) 
      throw std::exception();
    // that's to show it works
    if (blitz::kahan_sum(a) - blitz::kahan_sum(a.reverse(0)) > 4e-8) 
      throw std::exception();
  }


  ////////////////////////////////////////////////////////////////
  // Kahan sum functions
  {
    blitz::Array<real_t, 1> a(10), b(10);
    a = 2;
    b = a; //kahan_sum_5arg(a,a,a,a,a);
//    auto xpr = a+b;
//    auto itr = xpr.unwrap().begin();
  }
}

int main()
{
  main_tmpl<float>();
  main_tmpl<double>();
  main_tmpl<long double>();
}
