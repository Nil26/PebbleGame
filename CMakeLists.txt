cmake_minimum_required(VERSION 3.5)
project(PebbleGameTest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    SiteRP.cpp siteRP.h bond.h bond.cpp main.cpp)

add_executable(PebbleGameTest ${SOURCE_FILES})