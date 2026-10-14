#!/usr/bin/env bash
set -euo pipefail
: "${COPT_HOME:?Set COPT_HOME to your separately installed COPT 8.0.7 SDK}"
CXX="${CXX:-c++}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/build/mainline}"
mkdir -p "$BUILD_DIR"
FLAGS=(-std=c++20 -O3 -DNDEBUG -flto=1
  -DILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD=1
  -DILS_SP_FORWARD_ONLY_TWO_OPT=1
  -DILS_SP_TYPED_DELETE_ENTRY=0
  -DILS_SP_PACKED_DFS_MASK=1
  -I"$ROOT/include" -I"$ROOT/src"
  -isystem "$COPT_HOME/include" -isystem "$COPT_HOME/include/coptcpp_inc")
OBJECTS=()
for unit in core xml evaluator search pool copt_models solver; do
  "$CXX" "${FLAGS[@]}" -c "$ROOT/src/$unit.cpp" -o "$BUILD_DIR/$unit.o"
  OBJECTS+=("$BUILD_DIR/$unit.o")
done
for entry in main main_with_schedule; do
  "$CXX" "${FLAGS[@]}" -c "$ROOT/src/$entry.cpp" -o "$BUILD_DIR/$entry.o"
  output=ils-sp
  [[ "$entry" == main_with_schedule ]] && output=ils-sp-with-schedule
  "$CXX" -std=c++20 -O3 -DNDEBUG -flto=1 "${OBJECTS[@]}" "$BUILD_DIR/$entry.o" \
    -L"$COPT_HOME/lib" -Wl,-rpath,"$COPT_HOME/lib" \
    -lcopt_cpp -lcopt -l:libxxhash.so.0 -pthread -o "$BUILD_DIR/$output"
done
