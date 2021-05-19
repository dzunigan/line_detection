cmake_minimum_required(VERSION 3.0)

project(AG3line CXX)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")
# message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

## Set CXX optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wpedantic -Wall -Wextra")
# -mtune=native -march=native

## Set build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
endif(CMAKE_BUILD_TYPE STREQUAL "Release")
message(STATUS "Flags: " ${CMAKE_CXX_FLAGS})

## Dependencies
find_package(OpenCV REQUIRED)

include_directories(include ${PROJECT_SOURCE_DIR} ${OpenCV_INCLUDE_DIRS})

add_library(ag3line SHARED src/ag3line.cpp src/BaseLine.cpp src/fastatan2.cpp src/GeneralFuncs.cpp src/RegManager.cpp)
target_link_libraries(ag3line ${OpenCV_LIBS})
