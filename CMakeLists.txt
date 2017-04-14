cmake_minimum_required(VERSION 3.6)
project(proj_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES proj_2.c)

add_executable(proj_2 ${SOURCE_FILES})

find_library(M_LIB m X11)
find_library(CPGPLOT_LIB cpgplot pgplot png z)
target_link_libraries(proj_2 ${M_LIB})
target_link_libraries(proj_2 ${CPGPLOT_LIB})