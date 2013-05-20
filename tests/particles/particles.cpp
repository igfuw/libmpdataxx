#include <iostream>

#include "particles.hpp"

template <typename real_t, int thrust_device_system>
void particles<real_t, thrust_device_system>::func()
{
  std::cerr << "CUDA/OpenMP/CPP: " << thrust_device_system << std::endl;
}

template class particles<float, device_system_macro>;
