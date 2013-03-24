set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_package(CUDA QUIET)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) # -std=c++11 will not work

# trying to allow runtime choice betwee CUDA and OpenMP
set(tmp "")
set(tmp "${tmp};test_particles.cpp")
set(tmp "${tmp};particles_openmp.cpp")

if (CUDA_FOUND)
  set(tmp "${tmp};particles_cuda.cu")
  cuda_add_executable(test_particles ${tmp})
else()
  add_executable(test_particles ${tmp})
endif()
