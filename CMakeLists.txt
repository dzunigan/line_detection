cmake_minimum_required(VERSION 3.0)

project(line_detection CXX)
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

## Dependencies options
set(GFLAGS_PREFER_EXPORTED_GFLAGS_CMAKE_CONFIGURATION OFF CACHE BOOL "Use exported cmake config to find gflags?")

## Dependencies
find_package(Boost REQUIRED COMPONENTS filesystem system)
find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Gflags REQUIRED)
find_package(glog 0.5.0 REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(include ${PROJECT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GFLAGS_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})

## Header files
file(GLOB_RECURSE HEADER_FILES include/*.hpp include/*.h)
add_custom_target(header_files SOURCES ${HEADER_FILES})
message(STATUS "Header files: " ${HEADER_FILES})

add_executable(demo_tum src/demo_tum.cc src/ED.cc src/EDLines.cc src/NFA.cc src/FOVUndistorter.cc
  src/util/logging.cc src/util/math.cc src/util/string.cc src/util/timer.cc)
target_link_libraries(demo_tum ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${GLOG_LIBRARIES} ${OpenCV_LIBS})

add_executable(demo_dlsd src/demo_dlsd.cc src/ED.cc src/camera_models.cc src/polynomial.cc
  src/util/logging.cc src/util/math.cc src/util/random.cc src/util/string.cc src/util/timer.cc)
target_link_libraries(demo_dlsd glog::glog ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${OpenCV_LIBS})

add_executable(demo_ag3line src/demo_ag3line.cc src/FOVUndistorter.cc
  src/util/logging.cc src/util/math.cc src/util/random.cc src/util/string.cc src/util/timer.cc)
target_link_libraries(demo_ag3line glog::glog ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${OpenCV_LIBS} ${PROJECT_SOURCE_DIR}/thirdparty/AG3line/build/libag3line.so)

add_executable(demo_lsd src/demo_lsd.cc src/FOVUndistorter.cc
  src/util/logging.cc src/util/math.cc src/util/string.cc src/util/timer.cc)
target_link_libraries(demo_lsd ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${GLOG_LIBRARIES} ${OpenCV_LIBS} ${PROJECT_SOURCE_DIR}/thirdparty/LSD/build/liblsd.so)

add_executable(evaluate_tum src/evaluate_tum.cc src/dlsd.cc src/ED.cc src/EDLines.cc src/NFA.cc src/FOVUndistorter.cc src/camera_models.cc src/polynomial.cc
  src/util/logging.cc src/util/math.cc src/util/random.cc src/util/string.cc src/util/timer.cc)
target_link_libraries(evaluate_tum glog::glog ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${OpenCV_LIBS} ${PROJECT_SOURCE_DIR}/thirdparty/AG3line/build/libag3line.so ${PROJECT_SOURCE_DIR}/thirdparty/LSD/build/liblsd.so)

add_executable(rectify src/rectify.cc src/FOVUndistorter.cc src/util/logging.cc src/util/string.cc)
target_link_libraries(rectify glog::glog ${Boost_LIBRARIES} ${GFLAGS_LIBRARY} ${OpenCV_LIBS})
