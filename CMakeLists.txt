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

cmake_minimum_required(VERSION "4.3")

project(
  fcarouge-typed-linear-algebra
  DESCRIPTION "Typed Linear Algebra"
  HOMEPAGE_URL "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
  LANGUAGES "CXX"
  VERSION "0.4.0")

if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_CXX_STANDARD 26)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CPACK_PACKAGE_CHECKSUM "SHA512")
  set(CPACK_PACKAGE_CONTACT "francois.carouge@gmail.com")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
  set(CPACK_SOURCE_GENERATOR "TBZ2")

  include(CheckIPOSupported)
  include(CPack)
  include(CTest)
  include(FetchContent)
  include(GNUInstallDirs)

  include(support/support.cmake)

  add_subdirectory("benchmark")
  add_subdirectory("pkgconfig")
  add_subdirectory("sample")
  add_subdirectory("support")

  # Single-header amalgamations for tools that cannot resolve this project's
  # multi-file layout, such as Compiler Explorer's `#include`-by-URL feature.
  #
  # * NAME The generated header's file name, without extension.
  # * ENTRY The header to start amalgamating from.
  # * INCLUDE_DIRS Fallback roots to resolve quoted includes against, in
  #   order, after the including file's own directory.
  function(add_amalgamation NAME ENTRY)
    set(multiValueArgs INCLUDE_DIRS)
    cmake_parse_arguments(PARSE_ARGV 2 AMALGAMATION "" "" "${multiValueArgs}")

    set(OUTPUT_HEADER "${CMAKE_BINARY_DIR}/amalgamate/fcarouge/${NAME}.h")
    add_custom_command(
      OUTPUT "${OUTPUT_HEADER}"
      COMMAND
        "${CMAKE_COMMAND}" -D "INPUT=${ENTRY}" -D
        "INCLUDE_DIRS=${AMALGAMATION_INCLUDE_DIRS}" -D
        "OUTPUT=${OUTPUT_HEADER}" -P "${CMAKE_SOURCE_DIR}/support/amalgamate.cmake"
      DEPENDS "${ENTRY}"
      COMMENT "Generating amalgamated ${NAME}.h header"
      VERBATIM)
    add_custom_target("amalgamate_${NAME}" DEPENDS "${OUTPUT_HEADER}")
  endfunction(add_amalgamation)

  add_amalgamation(
    typed_linear_algebra
    "${CMAKE_SOURCE_DIR}/include/fcarouge/typed_linear_algebra.hpp"
    INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include")

  add_amalgamation(
    typed_linear_algebra_mp_units_eigen
    "${CMAKE_SOURCE_DIR}/support/mp_units_eigen/fcarouge/linalg.hpp"
    INCLUDE_DIRS
      "${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/support/eigen"
      "${CMAKE_SOURCE_DIR}/support/mp_units"
      "${CMAKE_SOURCE_DIR}/support/mp_units_eigen")

  add_custom_target(
    amalgamate DEPENDS amalgamate_typed_linear_algebra
                       amalgamate_typed_linear_algebra_mp_units_eigen)

  add_subdirectory("test")
endif()

add_subdirectory("cmake")
add_subdirectory("include")
