include(${CMAKE_SOURCE_DIR}/_EXTERNAL/MODULES/_customCMAKE.txt)


print("Checking compilers..." Boldmagenta)

SET(APPNAME ov)
# Project name
project(${APPNAME}_${VERSION} C CXX Fortran)

# Minimum version of CMake that can be used
cmake_minimum_required(VERSION 2.6)


print("" BoldYellow)
print("Finding packages..." Boldmagenta)



#SETTING BUILDING FLAGS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")





FILE(GLOB sources *.f *.c *.cpp)

add_library(${APPNAME} STATIC ${sources})

print("" BoldYellow)
print("Finishing..." Boldmagenta)
