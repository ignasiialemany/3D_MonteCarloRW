cmake_minimum_required(VERSION 3.5)
#set(CMAKE_C_COMPILER "/opt/homebrew/opt/gcc/bin/gcc-12")
#set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/gcc/bin/g++-12")
set(CMAKE_CXX_STANDARD 14)
project(3DRandomWalk)

set(CMAKE_BUILD_TYPE "Debug")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Choose your desired optimization level and floating-point settings
    set(OPTIMIZATION_LEVEL "-O3") 
    set(FP_CONSISTENCY_FLAGS "")

    # Apply the same optimization level and floating-point settings to Debug and Release configurations
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPTIMIZATION_LEVEL} ${FP_CONSISTENCY_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OPTIMIZATION_LEVEL} ${FP_CONSISTENCY_FLAGS}")
endif()

add_subdirectory(third_party/yaml-cpp)
add_subdirectory(third_party/Catch2-3.3.1)

# Set the name of the conda environment containing the required packages
set(CONDA_ENV_NAME "native")

# Set the paths to the conda environment
set(CONDA_PREFIX "/Users/ia4118/opt/anaconda3/envs/${CONDA_ENV_NAME}")
set(CONDA_INCLUDE_DIR "${CONDA_PREFIX}/include")
set(CONDA_LIBRARY_DIR "${CONDA_PREFIX}/lib")

#set(Eigen3_DIR "${CONDA_PREFIX}/share/eigen3/cmake")
set(CMAKE_PREFIX_PATH ${CONDA_PREFIX} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${PROJECT_SOURCE_DIR}/third_party/yaml-cpp")


find_package(Eigen3 3.3 REQUIRED)
message(STATUS "Eigen3 found at: ${Eigen3_DIR}")

find_package(yaml-cpp REQUIRED)

find_package(CGAL REQUIRED)
message(STATUS "cgal found at: ${CGAL_DIR}")

find_package(Boost REQUIRED COMPONENTS filesystem system)
message(STATUS "Boost found at: ${Boost_DIR}")

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found at: ${OpenMP_CXX_FLAGS}")

find_package(Catch2 REQUIRED)
message(STATUS "Catch2 found at: ${Catch2_DIR}")

find_package(matioCpp REQUIRED)
message(STATUS "matioCpp found at: ${matioCpp_DIR}")

set(GEOMETRY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Geometry")
set(SUBSTRATE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Substrate")
set(MONTECARLO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MonteCarlo")
set(MRI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MRI")
set(UTILITY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utility")

file(GLOB_RECURSE SOURCES "${GEOMETRY_DIR}/*.cpp" "${SUBSTRATE_DIR}/*.cpp" "${MONTECARLO_DIR}/*.cpp" "${MRI_DIR}/*.cpp" "${UTILITY_DIR}/*.cpp")
file(GLOB_RECURSE HEADERS "${GEOMETRY_DIR}/*.h" "${SUBSTRATE_DIR}/*.h" "${MONTECARLO_DIR}/*.h" "${MRI_DIR}/*.h" "${UTILITY_DIR}/*.h")

# Add the source files for your project
add_executable(myTests testing/oneGenerator.h testing/testing_simulation.cpp testing/testing_sequence.cpp testing/testing_polygon.cpp testing/testing_substrate.cpp ${SOURCES} ${HEADERS})
add_executable(3DRandomWalk run_sim.cpp ${SOURCES} ${HEADERS})

# Link against the yaml-cpp library
target_link_libraries(3DRandomWalk Eigen3::Eigen Catch2::Catch2WithMain OpenMP::OpenMP_CXX ${Boost_LIBRARIES} CGAL::CGAL matioCpp::matioCpp Boost::filesystem Boost::system yaml-cpp)
target_link_libraries(myTests Eigen3::Eigen Catch2::Catch2WithMain yaml-cpp OpenMP::OpenMP_CXX ${Boost_LIBRARIES} CGAL::CGAL matioCpp::matioCpp Boost::filesystem Boost::system)
