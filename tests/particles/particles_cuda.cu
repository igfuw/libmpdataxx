// including it first not to require pthread option to nvcc
//#include <blitz/array.h>

#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
#define device_system_macro cuda
#include "particles.cpp"
