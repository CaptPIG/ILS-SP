cmake_minimum_required(VERSION 3.24)

project(ILS_SP VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED COPT_HOME)
  if(DEFINED ENV{COPT_HOME})
    set(COPT_HOME "$ENV{COPT_HOME}")
  elseif(EXISTS "$ENV{HOME}/.local/opt/copt80")
    set(COPT_HOME "$ENV{HOME}/.local/opt/copt80")
  elseif(EXISTS "/opt/copt80")
    set(COPT_HOME "/opt/copt80")
  endif()
endif()

find_path(COPT_INCLUDE_DIR
  NAMES coptcpp.h
  HINTS "${COPT_HOME}/include"
  REQUIRED
)
find_library(COPT_CXX_LIBRARY
  NAMES copt_cpp
  HINTS "${COPT_HOME}/lib"
  REQUIRED
)
find_library(COPT_LIBRARY
  NAMES copt
  HINTS "${COPT_HOME}/lib"
  REQUIRED
)

find_package(Threads REQUIRED)
# Link a separately installed xxHash library; no bundled SDK or downloads.
find_library(ILS_SP_XXHASH_LIBRARY NAMES xxhash libxxhash.so.0 REQUIRED)

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

target_include_directories(ils_sp
  PUBLIC include
)
target_include_directories(ils_sp SYSTEM PUBLIC
  "${COPT_INCLUDE_DIR}" "${COPT_INCLUDE_DIR}/coptcpp_inc"
)
target_link_libraries(ils_sp
  PUBLIC "${COPT_CXX_LIBRARY}" "${COPT_LIBRARY}" Threads::Threads "${ILS_SP_XXHASH_LIBRARY}"
)
target_compile_options(ils_sp PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)

add_executable(ils-sp src/main.cpp)
target_link_libraries(ils-sp PRIVATE ils_sp)
