cmake_minimum_required(VERSION 3.24)

project(ILS_SP VERSION 1.0.0 LANGUAGES CXX)

option(ILS_SP_BUILD_TESTS "Build the C++ regression tests" ON)

if(NOT DEFINED GUROBI_HOME)
  if(DEFINED ENV{GUROBI_HOME})
    set(GUROBI_HOME "$ENV{GUROBI_HOME}")
  elseif(EXISTS "/opt/gurobi-current")
    set(GUROBI_HOME "/opt/gurobi-current")
  endif()
endif()

find_path(GUROBI_INCLUDE_DIR
  NAMES gurobi_c++.h
  HINTS "${GUROBI_HOME}/include"
  REQUIRED
)
find_library(GUROBI_CXX_LIBRARY
  NAMES gurobi_c++
  HINTS "${GUROBI_HOME}/lib"
  REQUIRED
)
find_library(GUROBI_LIBRARY
  NAMES gurobi130 gurobi120 gurobi110 gurobi100 gurobi
  HINTS "${GUROBI_HOME}/lib"
  REQUIRED
)

find_package(Threads REQUIRED)

add_library(ils_sp
  src/core.cpp
  src/xml.cpp
  src/evaluator.cpp
  src/search.cpp
  src/pool.cpp
  src/gurobi_models.cpp
  src/solver.cpp
)
add_library(ILS_SP::ils_sp ALIAS ils_sp)

target_compile_features(ils_sp PUBLIC cxx_std_20)
set_target_properties(ils_sp PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(ils_sp
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "${GUROBI_INCLUDE_DIR}"
)
target_link_libraries(ils_sp
  PUBLIC "${GUROBI_CXX_LIBRARY}" "${GUROBI_LIBRARY}" Threads::Threads
)
if(MSVC)
  target_compile_options(ils_sp PRIVATE /W4 /permissive-)
else()
  target_compile_options(
    ils_sp PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow
  )
endif()

add_executable(ils-sp src/main.cpp)
target_link_libraries(ils-sp PRIVATE ILS_SP::ils_sp)

if(ILS_SP_BUILD_TESTS)
  enable_testing()
  add_executable(ils-sp-tests tests/test_main.cpp)
  target_link_libraries(ils-sp-tests PRIVATE ILS_SP::ils_sp)
  add_test(NAME ils-sp-tests COMMAND ils-sp-tests "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(ils-sp-gurobi-tests tests/test_gurobi.cpp)
  target_link_libraries(ils-sp-gurobi-tests PRIVATE ILS_SP::ils_sp)
  add_test(NAME ils-sp-gurobi-tests COMMAND ils-sp-gurobi-tests)
endif()
