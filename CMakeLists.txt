cmake_minimum_required(VERSION 3.5)
set(CMAKE_C_COMPILER "/opt/homebrew/opt/gcc/bin/gcc-12")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/gcc/bin/g++-12")
set(CMAKE_CXX_STANDARD 14)
project(3DRandomWalk)

#set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_BUILD_TYPE Debug)
# Add compiler flags for debugging
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-function -Wno-error=deprecated-declarations -Wno-error=unused-variable -Wno-error=maybe-uninitialized -Wno-error=unused-parameter -Wno-error=unused-but-set-variable -Wno-error=unused-result -Wno-error=sign-compare -Wno-error=type-limits -Wno-error=parentheses -Wno-error=reorder -Wno-error=ignored-qualifiers -Wno-error=unused-local-typedefs -Wno-error=class-memaccess -Wno-error=nonnull -Wno-error=type-limits -Wno-error=misleading-indentation -std=c++17")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

# Add the yaml-cpp library
add_subdirectory(third_party/yaml-cpp)
add_subdirectory(third_party/Catch2-3.3.1)
add_subdirectory(third_party/matio-cpp-master)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem)
find_package(CGAL REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
include_directories(${CGAL_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(third_party/yaml-cpp/include)
include_directories(third_party/matio-cpp-master/include)

set(GEOMETRY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Geometry")
set(SUBSTRATE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Substrate")
set(MONTECARLO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MonteCarlo")
set(MRI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MRI")

file(GLOB_RECURSE SOURCES "${GEOMETRY_DIR}/*.cpp" "${SUBSTRATE_DIR}/*.cpp" "${MONTECARLO_DIR}/*.cpp" "${MRI_DIR}/*.cpp")
file(GLOB_RECURSE HEADERS "${GEOMETRY_DIR}/*.h" "${SUBSTRATE_DIR}/*.h" "${MONTECARLO_DIR}/*.h" "${MRI_DIR}/*.h")

# Add the source files for your project
add_executable(myTests testing/testing_polygon.cpp testing/testing_substrate.cpp ${SOURCES} ${HEADERS})
add_executable(3DRandomWalk run_sim.cpp ${SOURCES} ${HEADERS})

# Link against the yaml-cpp library
target_link_libraries(3DRandomWalk Eigen3::Eigen yaml-cpp OpenMP::OpenMP_CXX ${Boost_LIBRARIES} CGAL::CGAL matioCpp::matioCpp)
target_link_libraries(myTests Eigen3::Eigen Catch2::Catch2WithMain yaml-cpp OpenMP::OpenMP_CXX ${Boost_LIBRARIES} CGAL::CGAL matioCpp::matioCpp Boost::filesystem Boost::system)



