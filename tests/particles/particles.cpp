#include <iostream>

#include "particles.hpp"

#include <thrust/device_vector.h>

template <typename real_t, int thrust_device_system>
void particles<real_t, thrust_device_system>::func()
{
#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP) 
  assert(device_system_macro == openmp);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA) 
  assert(device_system_macro == cuda);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP) 
  assert(device_system_macro == cpp);
#endif

  std::cerr << "CUDA/OpenMP/CPP: " << thrust_device_system << std::endl;
  thrust::device_vector<real_t> vec(1024*1024);
}

template class particles<float, device_system_macro>;
