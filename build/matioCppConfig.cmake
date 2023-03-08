set(matioCpp_VERSION 0.2.2)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was matioCppConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################

#### Expanded from @PACKAGE_DEPENDENCIES@ by install_basic_package_files() ####

include(CMakeFindDependencyMacro)
set(CMAKE_MODULE_PATH_BK_matioCpp ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${PACKAGE_PREFIX_DIR}/share/matioCpp/cmake)
find_package(MATIO)
find_package(visit_struct)
find_package(Eigen3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_BK_matioCpp})

###############################################################################


if(NOT TARGET matioCpp::matioCpp)
  include("${CMAKE_CURRENT_LIST_DIR}/matioCppTargets.cmake")
endif()

# Compatibility
get_property(matioCpp_matioCpp_INCLUDE_DIR TARGET matioCpp::matioCpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

set(matioCpp_LIBRARIES matioCpp::matioCpp)
set(matioCpp_INCLUDE_DIRS "${matioCpp_matioCpp_INCLUDE_DIR}")
list(REMOVE_DUPLICATES matioCpp_INCLUDE_DIRS)


