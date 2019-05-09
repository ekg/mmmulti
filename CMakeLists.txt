# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(mmmultimap)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

# Use openmp for parallelism, but it's configured differently on OSX
find_package(OpenMP)
if (OPENMP_FOUND)
  if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
  endif()
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# TODO: We're using INSTALL_DIR very wrong. We *should* be actually installing
# the external projects into their prefixes and working with the installed
# files. Instead we're building but not installing them and trying to work with
# the non-installed build trees.
# 
# Hence the blanked out INSTALL_COMMANDs to suppress the install step.
#
# We need to NOT blank out UPDATE_COMMAND or we can never change the Git revision we point to.
# The cost of this is that we have to re-configure on every build if we do update.

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "d52aa9a71513d132e30c09491b5899af449ebb94"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  GIT_REPOSITORY "https://github.com/vgteam/DYNAMIC.git"
  GIT_TAG "615d8be5276bcd4c5a3d8e31679b4f8e81b2eefc"
  # we don't actually install dynamic... it's header only
  #CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(dynamic INSTALL_DIR)
set(dynamic_INCLUDE "${INSTALL_DIR}/src/dynamic/include")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  GIT_REPOSITORY "https://github.com/SaschaWitt/ips4o.git"
  GIT_TAG "bff3ccf0bf349497f2bb10f825d160b792236367"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

set(CMAKE_BUILD_TYPE Release)

# set up our target executable and specify its dependencies and includes
add_executable(mmmultimap
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  ${CMAKE_SOURCE_DIR}/src/buffer.cpp
  ${CMAKE_SOURCE_DIR}/src/threads.cpp
  )
add_dependencies(mmmultimap sdsl-lite)
add_dependencies(mmmultimap dynamic)
add_dependencies(mmmultimap ips4o)
add_dependencies(mmmultimap tayweeargs)
target_include_directories(mmmultimap PUBLIC
  "${CMAKE_SOURCE_DIR}/src"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${dynamic_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${tayweeargs_INCLUDE}")
target_link_libraries(mmmultimap
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "-latomic"
#  "-ltcmalloc_minimal"
  )

#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltcmalloc")

if (APPLE)
elseif (TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
