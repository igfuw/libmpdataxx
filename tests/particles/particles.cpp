#include <iostream>

#include "particles.hpp"

#include <thrust/device_vector.h>

template <typename real_t, int thrust_device_system>
void particles<real_t, thrust_device_system>::func()
{
  std::cerr << "CUDA/OpenMP/CPP: " << thrust_device_system << std::endl;
  thrust::device_vector<real_t> vec(1024*1024);
}

template class particles<float, device_system_macro>;
