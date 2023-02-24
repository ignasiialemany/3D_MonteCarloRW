cmake_minimum_required(VERSION 3.5)
project(3DRandomWalk)

# Set the C++ standard to C++11
set(CMAKE_CXX_STANDARD 17)

# Add the yaml-cpp library
add_subdirectory(third_party/yaml-cpp)
add_subdirectory(third_party/Catch2-3.3.1)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Set the include directories
include_directories(third_party/yaml-cpp/include)
include_directories(third_party/eigen-3.4.0)

# Add the source files for your project
add_executable(myTests testing/testing_walkers.cpp testing/testing_simulation.cpp MonteCarlo/walkers.cpp MonteCarlo/walkers.h MonteCarlo/simulation.cpp MonteCarlo/simulation.h)
add_executable(3DRandomWalk run_sim.cpp MonteCarlo/walkers.cpp MonteCarlo/walkers.h MonteCarlo/simulation.cpp MonteCarlo/simulation.h)

# Link against the yaml-cpp library
target_link_libraries(3DRandomWalk Eigen3::Eigen yaml-cpp)
target_link_libraries(myTests PRIVATE Catch2::Catch2WithMain)
