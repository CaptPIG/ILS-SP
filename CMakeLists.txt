cmake_minimum_required(VERSION 3.24)

project(ILS_SP VERSION 1.0.0 LANGUAGES CXX)

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(FATAL_ERROR "This frozen mainline requires Linux filesystem APIs.")
endif()

if(NOT DEFINED COPT_HOME AND DEFINED ENV{COPT_HOME})
  set(COPT_HOME "$ENV{COPT_HOME}")
endif()
find_path(COPT_INCLUDE_DIR NAMES coptcpp.h HINTS "${COPT_HOME}/include" REQUIRED)
find_library(COPT_CXX_LIBRARY NAMES copt_cpp HINTS "${COPT_HOME}/lib" REQUIRED)
find_library(COPT_LIBRARY NAMES copt HINTS "${COPT_HOME}/lib" REQUIRED)
find_library(ILS_SP_XXHASH_LIBRARY NAMES xxhash libxxhash.so.0 REQUIRED)
find_package(Threads REQUIRED)

add_library(ils_sp
  src/core.cpp
  src/xml.cpp
  src/evaluator.cpp
  src/search.cpp
  src/pool.cpp
  src/copt_models.cpp
  src/solver.cpp
)
add_library(ILS_SP::ils_sp ALIAS ils_sp)
target_compile_features(ils_sp PUBLIC cxx_std_20)
set_target_properties(ils_sp PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(ils_sp PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_include_directories(ils_sp SYSTEM PUBLIC
  "${COPT_INCLUDE_DIR}" "${COPT_INCLUDE_DIR}/coptcpp_inc"
)
target_compile_definitions(ils_sp PUBLIC
  ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD=1
  ILS_SP_FORWARD_ONLY_TWO_OPT=1
  ILS_SP_TYPED_DELETE_ENTRY=0
  ILS_SP_PACKED_DFS_MASK=1
)
target_link_libraries(ils_sp PUBLIC
  "${COPT_CXX_LIBRARY}" "${COPT_LIBRARY}" Threads::Threads "${ILS_SP_XXHASH_LIBRARY}"
)
target_compile_options(ils_sp PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)

add_executable(ils-sp src/main.cpp)
target_link_libraries(ils-sp PRIVATE ILS_SP::ils_sp)
add_executable(ils-sp-with-schedule src/main_with_schedule.cpp)
target_link_libraries(ils-sp-with-schedule PRIVATE ILS_SP::ils_sp)

# The frozen R4 recipe used GCC 14.2.0, -O3 -DNDEBUG -flto=1.
# Keep LTO optional for other Linux toolchains and development builds.
option(ILS_SP_ENABLE_LTO "Enable link-time optimization for optimized builds" OFF)
if(ILS_SP_ENABLE_LTO)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(ils_sp PUBLIC -flto=1)
    target_link_options(ils_sp PUBLIC -flto=1)
  else()
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
    if(NOT ipo_supported)
      message(FATAL_ERROR "Link-time optimization is unavailable: ${ipo_error}")
    endif()
    set_property(TARGET ils_sp ils-sp ils-sp-with-schedule PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endif()

# Existing Gurobi implementation/test files remain in the repository for
# historical context; this source-only mainline build does not compile them.
