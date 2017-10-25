cmake_minimum_required(VERSION 3.5)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(PebbleGameTest)

#option(Testing "Enable testing" ON)
#include(EnableGtests) #defined in ./cmake

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -v")

set(SOURCE_FILES
    SiteRP.cpp bond.cpp main.cpp SiteRP_cont.cpp )

#include_directories(${GTEST_INCLUDE_DIR})
add_executable(PebbleGameTest ${SOURCE_FILES})

#set(unittest_src
#   SiteRP_cont_test
#        )

#enable_testing()
#foreach(test ${unittest_src})
#    add_gtest(${test})
#endforeach(test)
