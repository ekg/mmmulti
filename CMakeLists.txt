# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(mmmulti)
# We build using c++14
set(CMAKE_CXX_STANDARD 17)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcx16 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcx16 -g")
# for debugging, unset and run with ASAN_OPTIONS=detect_leaks=1:symbolize=1 LSAN_OPTIONS=verbosity=1:log_threads=1
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -mcx16 -g -fsanitize=address")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -mcx16 -g -fsanitize=address")

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
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sdsl-lite"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/DYNAMIC"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(dynamic SOURCE_DIR)
set(dynamic_INCLUDE "${SOURCE_DIR}/include")

ExternalProject_Add(hopscotch_map
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/hopscotch-map"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(hopscotch_map SOURCE_DIR)
set(hopscotch_map_INCLUDE "${SOURCE_DIR}/include")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
add_subdirectory(${CMAKE_SOURCE_DIR}/deps/ips4o)

# atomic queue
ExternalProject_Add(atomicqueue
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomic_queue/include/atomic_queue"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicqueue SOURCE_DIR)
set(atomicqueue_INCLUDE "${SOURCE_DIR}/")

# paryfor parallel_for
ExternalProject_Add(paryfor
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/paryfor"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(paryfor SOURCE_DIR)
set(paryfor_INCLUDE "${SOURCE_DIR}")

# taywee's C++ args library, header only
set(tayweeargs_INCLUDE deps/args)

add_subdirectory(deps/mio)

#set(CMAKE_BUILD_TYPE Release)

# set up our target executable and specify its dependencies and includes
add_executable(mmmulti
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  )
add_dependencies(mmmulti sdsl-lite)
add_dependencies(mmmulti dynamic)
add_dependencies(mmmulti atomicqueue)
add_dependencies(mmmulti paryfor)

target_include_directories(mmmulti PUBLIC
  "${CMAKE_SOURCE_DIR}/src"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${dynamic_INCLUDE}"
  "${atomicqueue_INCLUDE}"
  "${paryfor_INCLUDE}"
  "${tayweeargs_INCLUDE}")
  
# macOS doesn't want you to link in libatomic this way
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_link_libraries(mmmulti
    "${sdsl-lite_LIB}/libsdsl.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
    Threads::Threads
    ips4o
    mio
    )
elseif (TRUE)
  target_link_libraries(mmmulti
    "${sdsl-lite_LIB}/libsdsl.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
    "-latomic"
    Threads::Threads
    ips4o
    mio
    )
endif()


if (APPLE)
elseif (TRUE)
  # this was hard to track down
  # https://stackoverflow.com/questions/35116327/when-g-static-link-pthread-cause-segmentation-fault-why
  #set(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
endif()



