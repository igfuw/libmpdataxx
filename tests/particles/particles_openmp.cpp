#if defined(_OPENMP) 
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
#  define device_system_macro openmp
#  include "particles.cpp"
#endif
// TODO else ..._SYSTEM_CPP when it will finally be supported
