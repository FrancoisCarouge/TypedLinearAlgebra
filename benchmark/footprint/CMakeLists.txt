#[[ Typed Linear Algebra
Version 0.4.0
https://github.com/FrancoisCarouge/TypedLinearAlgebra

SPDX-License-Identifier: Unlicense

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org> ]]

if(NOT BUILD_TESTING)
  return()
endif()

# Add a given embedded footprint measurement.
#
# Cross-compiled, the test passes when the operation's step entry point, linked
# alone, measures the expected sizes in bytes. Without a variant, the typed
# matrix and its freestanding storage are measured against the same expected
# sizes: the typing costs zero bytes. Otherwise, the configuration is only
# compiled. The
# configuration `<rows>x<columns>` is configured from the `configuration.cpp`
# template. The tests are named after the directory they are declared in, their
# configuration, operation, optimization, and variant: `typed` or `storage`.
#
# * NAME The configuration: `<rows>x<columns>`.
# * OPERATION The operation of the entry points.
# * OPTIMIZATION The optimization level: `Os` or `O2`.
# * VARIANT The only variant measured: `typed` or `storage`, for an operation
#   without a storage counterpart, or whose typing costs bytes.
# * FLASH The expected read-only memory size: `.text`, `.rodata`, `.data`.
# * RAM The expected read-write memory size: `.data`, `.bss`.
# * STACK The expected worst-case stack size.
function(footprint NAME OPERATION OPTIMIZATION)
  cmake_parse_arguments(PARSE_ARGV 3 EXPECTED "" "VARIANT;FLASH;RAM;STACK" "")

  if(NOT NAME MATCHES "^([0-9]+)x([0-9]+)$")
    message(FATAL_ERROR "Unknown footprint configuration: ${NAME}")
  endif()
  set(ROWS ${CMAKE_MATCH_1})
  set(COLUMNS ${CMAKE_MATCH_2})
  configure_file("configuration.cpp" "${NAME}.cpp")

  if(NOT CMAKE_CROSSCOMPILING)
    if(NOT TARGET
       typed_linear_algebra_freestanding_footprint_${NAME}_build_driver)
      build("${NAME}" BACKENDS "freestanding")
      # The entry points' `[[gnu::flatten]]` shapes the cross-compiled image
      # only: MSVC does not recognize it (C5030).
      target_compile_options(
        typed_linear_algebra_freestanding_footprint_${NAME}_build_driver
        PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/wd5030>)
    endif()
    return()
  endif()

  set(MEMORY typed_linear_algebra_footprint_memory_${OPTIMIZATION})
  if(NOT TARGET ${MEMORY})
    add_library(${MEMORY} OBJECT "memory.cpp")
    target_compile_options(
      ${MEMORY}
      PRIVATE "-${OPTIMIZATION}" "-fno-tree-loop-distribute-patterns")
    target_link_libraries(
      ${MEMORY}
      PRIVATE typed_linear_algebra_footprint_options)
  endif()

  set(CONFIGURATION typed_linear_algebra_footprint_${NAME}_${OPTIMIZATION})
  if(NOT TARGET ${CONFIGURATION})
    add_library(${CONFIGURATION} OBJECT "${NAME}.cpp")
    target_compile_options(${CONFIGURATION} PRIVATE "-${OPTIMIZATION}")
    target_link_libraries(
      ${CONFIGURATION}
      PRIVATE
        tlinalg typed_linear_algebra_freestanding
        typed_linear_algebra_footprint_options)
  endif()

  set(VARIANTS "typed" "storage")
  if(DEFINED EXPECTED_VARIANT)
    if(NOT EXPECTED_VARIANT IN_LIST VARIANTS)
      message(FATAL_ERROR "Unknown footprint variant: ${EXPECTED_VARIANT}")
    endif()
    set(VARIANTS ${EXPECTED_VARIANT})
  endif()

  foreach(VARIANT IN LISTS VARIANTS)
    string(
      TOLOWER "${NAME}_${OPERATION}_${OPTIMIZATION}_${VARIANT}" MEASUREMENT)
    set(IMAGE typed_linear_algebra_footprint_${MEASUREMENT})
    set(ROW "${CMAKE_CURRENT_BINARY_DIR}/row/${MEASUREMENT}.html")
    add_executable(${IMAGE} $<TARGET_OBJECTS:${CONFIGURATION}>
                   $<TARGET_OBJECTS:${MEMORY}>)
    set_target_properties(${IMAGE} PROPERTIES LINKER_LANGUAGE "CXX")
    target_link_options(
      ${IMAGE}
      PRIVATE
        "-nostdlib" "LINKER:--gc-sections"
        "LINKER:--entry=${VARIANT}_${OPERATION}_step"
        "LINKER:--defsym=memcpy=footprint_copy"
        "LINKER:--defsym=memset=footprint_fill")
    add_test(
      NAME typed_linear_algebra_freestanding_footprint_${MEASUREMENT}_pass
      COMMAND
        ${CMAKE_COMMAND}
        "-DNM=${CMAKE_NM}"
        "-DREADELF=${CMAKE_READELF}"
        "-DIMAGE=$<TARGET_FILE:${IMAGE}>"
        "-DENTRY_USAGE=$<TARGET_OBJECTS:${CONFIGURATION}>"
        "-DLEAF_USAGE=$<TARGET_OBJECTS:${MEMORY}>"
        "-DNAME=${NAME}"
        "-DPOLICY=${VARIANT}_${OPERATION}"
        "-DOPTIMIZATION=${OPTIMIZATION}"
        "-DFLASH=${EXPECTED_FLASH}"
        "-DRAM=${EXPECTED_RAM}"
        "-DSTACK=${EXPECTED_STACK}"
        "-DROW=${ROW}"
        -P
        "${CMAKE_CURRENT_SOURCE_DIR}/footprint.cmake")
    set_tests_properties(
      typed_linear_algebra_freestanding_footprint_${MEASUREMENT}_pass
      PROPERTIES FIXTURES_REQUIRED "typed_linear_algebra_footprint_page")
    set_property(DIRECTORY APPEND PROPERTY FOOTPRINT_ROWS "${ROW}")
  endforeach()
endfunction(footprint)

# Add the documentation page test, assembling the rows of every declared
# measurement after them, passed or failed.
function(footprint_page)
  get_property(ROWS DIRECTORY PROPERTY FOOTPRINT_ROWS)
  add_test(
    NAME typed_linear_algebra_freestanding_footprint_page_pass
    COMMAND
      ${CMAKE_COMMAND}
      "-DROWS=${ROWS}"
      "-DTEMPLATE=${PROJECT_SOURCE_DIR}/documentation/html/footprint/index.html.in"
      "-DPAGE=${CMAKE_BINARY_DIR}/footprint/index.html"
      "-DCOMPILER=${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
      "-DFLAGS=${CMAKE_CXX_FLAGS}"
      "-DVERSION=${PROJECT_VERSION}"
      -P
      "${CMAKE_CURRENT_SOURCE_DIR}/footprint_page.cmake")
  set_tests_properties(
    typed_linear_algebra_freestanding_footprint_page_pass
    PROPERTIES FIXTURES_CLEANUP "typed_linear_algebra_footprint_page")
endfunction(footprint_page)

if(CMAKE_CROSSCOMPILING)
  add_library(typed_linear_algebra_footprint_options INTERFACE)
  target_compile_options(
    typed_linear_algebra_footprint_options
    INTERFACE "-ffunction-sections" "-fdata-sections" "-fstack-usage")
  cmake_language(DEFER CALL footprint_page)
else()
  add_library(typed_linear_algebra_footprint_memory OBJECT "memory.cpp")
  set_target_properties(
    typed_linear_algebra_footprint_memory
    PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

# Changing a criterion requires its reason in the commit message. The typed
# product copies its result once more than its storage, except for a 1x1
# matrix, multiplied as a scalar: its variants are measured apart.
footprint("1x1" "addition" "O2" FLASH 24 RAM 4 STACK 0)
footprint("1x1" "addition" "Os" FLASH 28 RAM 4 STACK 0)
footprint("1x1" "equal_to" "O2" VARIANT "typed" FLASH 24 RAM 0 STACK 0)
footprint("1x1" "equal_to" "Os" VARIANT "typed" FLASH 24 RAM 0 STACK 0)
footprint("1x1" "left_scale" "O2" FLASH 20 RAM 4 STACK 0)
footprint("1x1" "left_scale" "Os" FLASH 24 RAM 4 STACK 0)
footprint("1x1" "negation" "O2" FLASH 24 RAM 4 STACK 0)
footprint("1x1" "negation" "Os" FLASH 24 RAM 4 STACK 0)
footprint("1x1" "product" "O2" VARIANT "storage" FLASH 32 RAM 4 STACK 0)
footprint("1x1" "product" "O2" VARIANT "typed" FLASH 24 RAM 4 STACK 0)
footprint("1x1" "product" "Os" VARIANT "storage" FLASH 32 RAM 4 STACK 0)
footprint("1x1" "product" "Os" VARIANT "typed" FLASH 28 RAM 4 STACK 0)
footprint("1x1" "scale" "O2" FLASH 20 RAM 4 STACK 0)
footprint("1x1" "scale" "Os" FLASH 24 RAM 4 STACK 0)
footprint("1x1" "subtraction" "O2" FLASH 24 RAM 4 STACK 0)
footprint("1x1" "subtraction" "Os" FLASH 28 RAM 4 STACK 0)
footprint("2x3" "addition" "O2" FLASH 108 RAM 24 STACK 0)
footprint("2x3" "addition" "Os" FLASH 108 RAM 24 STACK 0)
footprint("2x3" "equal_to" "O2" VARIANT "typed" FLASH 148 RAM 0 STACK 0)
footprint("2x3" "equal_to" "Os" VARIANT "typed" FLASH 144 RAM 0 STACK 0)
footprint("2x3" "left_scale" "O2" FLASH 84 RAM 24 STACK 0)
footprint("2x3" "left_scale" "Os" FLASH 84 RAM 24 STACK 0)
footprint("2x3" "negation" "O2" FLASH 84 RAM 24 STACK 0)
footprint("2x3" "negation" "Os" FLASH 84 RAM 24 STACK 0)
footprint("2x3" "product" "O2" VARIANT "storage" FLASH 164 RAM 16 STACK 40)
footprint("2x3" "product" "O2" VARIANT "typed" FLASH 192 RAM 16 STACK 80)
footprint("2x3" "product" "Os" VARIANT "storage" FLASH 164 RAM 16 STACK 48)
footprint("2x3" "product" "Os" VARIANT "typed" FLASH 188 RAM 16 STACK 80)
footprint("2x3" "scale" "O2" FLASH 84 RAM 24 STACK 0)
footprint("2x3" "scale" "Os" FLASH 84 RAM 24 STACK 0)
footprint("2x3" "subtraction" "O2" FLASH 108 RAM 24 STACK 0)
footprint("2x3" "subtraction" "Os" FLASH 108 RAM 24 STACK 0)
footprint("3x3" "addition" "O2" FLASH 152 RAM 36 STACK 0)
footprint("3x3" "addition" "Os" FLASH 156 RAM 36 STACK 0)
footprint("3x3" "equal_to" "O2" VARIANT "typed" FLASH 220 RAM 0 STACK 0)
footprint("3x3" "equal_to" "Os" VARIANT "typed" FLASH 216 RAM 0 STACK 0)
footprint("3x3" "left_scale" "O2" FLASH 120 RAM 36 STACK 0)
footprint("3x3" "left_scale" "Os" FLASH 120 RAM 36 STACK 0)
footprint("3x3" "negation" "O2" FLASH 120 RAM 36 STACK 0)
footprint("3x3" "negation" "Os" FLASH 120 RAM 36 STACK 0)
footprint("3x3" "product" "O2" VARIANT "storage" FLASH 144 RAM 36 STACK 64)
footprint("3x3" "product" "O2" VARIANT "typed" FLASH 172 RAM 36 STACK 96)
footprint("3x3" "product" "Os" VARIANT "storage" FLASH 140 RAM 36 STACK 64)
footprint("3x3" "product" "Os" VARIANT "typed" FLASH 156 RAM 36 STACK 96)
footprint("3x3" "scale" "O2" FLASH 120 RAM 36 STACK 0)
footprint("3x3" "scale" "Os" FLASH 120 RAM 36 STACK 0)
footprint("3x3" "subtraction" "O2" FLASH 152 RAM 36 STACK 0)
footprint("3x3" "subtraction" "Os" FLASH 156 RAM 36 STACK 0)
footprint("4x1" "addition" "O2" FLASH 76 RAM 16 STACK 0)
footprint("4x1" "addition" "Os" FLASH 76 RAM 16 STACK 0)
footprint("4x1" "equal_to" "O2" VARIANT "typed" FLASH 96 RAM 0 STACK 0)
footprint("4x1" "equal_to" "Os" VARIANT "typed" FLASH 96 RAM 0 STACK 0)
footprint("4x1" "left_scale" "O2" FLASH 60 RAM 16 STACK 0)
footprint("4x1" "left_scale" "Os" FLASH 60 RAM 16 STACK 0)
footprint("4x1" "negation" "O2" FLASH 60 RAM 16 STACK 0)
footprint("4x1" "negation" "Os" FLASH 60 RAM 16 STACK 0)
footprint("4x1" "product" "O2" VARIANT "storage" FLASH 128 RAM 64 STACK 80)
footprint("4x1" "product" "O2" VARIANT "typed" FLASH 172 RAM 64 STACK 144)
footprint("4x1" "product" "Os" VARIANT "storage" FLASH 120 RAM 64 STACK 80)
footprint("4x1" "product" "Os" VARIANT "typed" FLASH 144 RAM 64 STACK 144)
footprint("4x1" "scale" "O2" FLASH 60 RAM 16 STACK 0)
footprint("4x1" "scale" "Os" FLASH 60 RAM 16 STACK 0)
footprint("4x1" "subtraction" "O2" FLASH 76 RAM 16 STACK 0)
footprint("4x1" "subtraction" "Os" FLASH 76 RAM 16 STACK 0)
footprint("4x4" "addition" "O2" FLASH 276 RAM 64 STACK 8)
footprint("4x4" "addition" "Os" FLASH 276 RAM 64 STACK 8)
footprint("4x4" "equal_to" "O2" VARIANT "typed" FLASH 388 RAM 0 STACK 0)
footprint("4x4" "equal_to" "Os" VARIANT "typed" FLASH 384 RAM 0 STACK 0)
footprint("4x4" "left_scale" "O2" FLASH 212 RAM 64 STACK 8)
footprint("4x4" "left_scale" "Os" FLASH 212 RAM 64 STACK 8)
footprint("4x4" "negation" "O2" FLASH 204 RAM 64 STACK 0)
footprint("4x4" "negation" "Os" FLASH 204 RAM 64 STACK 0)
footprint("4x4" "product" "O2" VARIANT "storage" FLASH 160 RAM 64 STACK 88)
footprint("4x4" "product" "O2" VARIANT "typed" FLASH 200 RAM 64 STACK 152)
footprint("4x4" "product" "Os" VARIANT "storage" FLASH 144 RAM 64 STACK 88)
footprint("4x4" "product" "Os" VARIANT "typed" FLASH 172 RAM 64 STACK 152)
footprint("4x4" "scale" "O2" FLASH 212 RAM 64 STACK 8)
footprint("4x4" "scale" "Os" FLASH 212 RAM 64 STACK 8)
footprint("4x4" "subtraction" "O2" FLASH 276 RAM 64 STACK 8)
footprint("4x4" "subtraction" "Os" FLASH 276 RAM 64 STACK 8)
