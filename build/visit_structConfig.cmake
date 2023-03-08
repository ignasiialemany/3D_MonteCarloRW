set(visit_struct_VERSION 0.0.1)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was visit_structConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################



if(NOT TARGET visit_struct::visit_struct)
  include("${CMAKE_CURRENT_LIST_DIR}/visit_structTargets.cmake")
endif()

# Compatibility
get_property(visit_struct_visit_struct_INCLUDE_DIR TARGET visit_struct::visit_struct PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

set(visit_struct_LIBRARIES visit_struct::visit_struct)
set(visit_struct_INCLUDE_DIRS "${visit_struct_visit_struct_INCLUDE_DIR}")
list(REMOVE_DUPLICATES visit_struct_INCLUDE_DIRS)


