project(MIA_2012_Vessel)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${OT_CMAKE_MODULE_PATH} ${OT_ICG_CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

make_directory(${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

if (WIN32)
  # for msvc 8 deprecated warning removal
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
  add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE_GLOBALS)
  add_definitions(-DWIN32 -D_WIN32)
  add_definitions(-DNOMINMAX)
  set(CMAKE_DEBUG_POSTFIX "D")
endif(WIN32)

# ITK:
find_package(ITK REQUIRED)
if(UNIX)
  set(ITK_REQUIRED_CXX_FLAGS "-Wall")
endif()
include(${ITK_USE_FILE})

# for openmp
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -g")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

#SET(CMAKE_BUILD_TYPE "Release")


include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/source)

add_executable(miaKUSheet1 source/main.cpp source/VolumeHandler.h source/VolumeHandler.cpp source/Definitions.h source/TVL1.cpp source/LungSegmentation.cpp)

target_link_libraries(miaKUSheet1 ${ITK_LIBRARIES})

