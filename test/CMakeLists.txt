cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
enable_language(CXX)
project (gf2sim)

add_definitions(-std=c++11)

#set(CMAKE_CXX_COMPILER /opt/cray/craype/2.4.2/bin/CC)
#set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options system serialization)
find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params gf)
find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)
option(Testing "Enable testing" ON)
include(EnableGtests) #defined in ./cmake

message("alps libraries" ${ALPSCore_LIBRARIES})
message("boost libraries" ${Boost_LIBRARIES})

include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${NFFT3_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS ")

add_library(bdmc
bdmc_hash.cpp
bdmc_graph.cpp
bdmc_vertex.cpp
bdmc_dot.cpp
bdmc_pseudo.cpp
)

target_link_libraries(bdmc
${ALPSCore_LIBRARIES}
${Boost_LIBRARIES}
${MPI_CXX_LIBRARIES}
${GSL_LIBRARIES}
${FFTW3_LIBRARIES}
${LAPACK_LIBRARIES}
)

add_executable(bdmcsim bdmc_main.cpp)
target_link_libraries(bdmcsim bdmc)

include_directories(${GTEST_INCLUDE_DIR})
set(unittest_src
bdmc_dot_test
bdmc_pseudo_test
bdmc_graph_test
bdmc_hash_test
bdmc_vertex_test)

enable_testing()
list(APPEND LINK_ALL bdmc)
foreach(test ${unittest_src})
    add_gtest(${test})
endforeach(test)
