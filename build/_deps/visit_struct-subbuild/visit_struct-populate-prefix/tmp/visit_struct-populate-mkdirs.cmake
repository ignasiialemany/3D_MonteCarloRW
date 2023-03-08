# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-src"
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-build"
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix"
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/tmp"
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp"
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src"
  "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/ia4118/CLionProjects/3D_RandomWalk/build/_deps/visit_struct-subbuild/visit_struct-populate-prefix/src/visit_struct-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
