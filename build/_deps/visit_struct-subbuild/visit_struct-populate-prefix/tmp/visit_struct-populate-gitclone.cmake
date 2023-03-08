# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitclone-lastrun.txt" AND EXISTS "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitinfo.txt" AND
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitclone-lastrun.txt" IS_NEWER_THAN "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git" 
            clone --no-checkout --config "advice.detachedHead=false" "https://github.com/ami-iit/visit_struct" "visit_struct-src"
    WORKING_DIRECTORY "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/ami-iit/visit_struct'")
endif()

execute_process(
  COMMAND "/usr/bin/git" 
          checkout "47bc6a3aa7588a1f4db39579a0b6812569a76b56" --
  WORKING_DIRECTORY "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '47bc6a3aa7588a1f4db39579a0b6812569a76b56'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-src"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitinfo.txt" "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/visit_struct-populate-gitclone-lastrun.txt'")
endif()
