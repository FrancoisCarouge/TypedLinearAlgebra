#[[ Typed Linear Algebra
Version 0.1.0
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

FetchContent_Declare(
  mdspan
  GIT_REPOSITORY "https://github.com/kokkos/mdspan.git"
  GIT_TAG "stable"
  FIND_PACKAGE_ARGS NAMES mdspan)
FetchContent_MakeAvailable(mdspan)

add_library(std::mdspan ALIAS mdspan)

FetchContent_Declare(
  linalg
  GIT_REPOSITORY "https://github.com/kokkos/stdBLAS.git"
  GIT_TAG "main"
  FIND_PACKAGE_ARGS NAMES linalg)
FetchContent_MakeAvailable(linalg)

add_library(typed_linear_algebra_kokkos INTERFACE)
target_sources(
  typed_linear_algebra_kokkos
  INTERFACE FILE_SET
            "typed_linear_algebra_headers"
            TYPE
            "HEADERS"
            FILES
            "fcarouge/linalg"
            "fcarouge/mdspan")
target_include_directories(typed_linear_algebra_kokkos INTERFACE "fcarouge")
target_link_libraries(typed_linear_algebra_kokkos INTERFACE linalg)

target_include_directories(
  typed_linear_algebra_kokkos
  INTERFACE ${linalg_BINARY_DIR}/include/experimental
            ${linalg_SOURCE_DIR}/include/experimental)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(
    typed_linear_algebra_kokkos INTERFACE "-Wno-unused-local-typedefs"
                                          "-Wno-shadow" "-Wno-unused-parameter")
endif()

target_compile_definitions(
  typed_linear_algebra_kokkos
  INTERFACE "MDSPAN_IMPL_STANDARD_NAMESPACE=Kokkos"
            "MDSPAN_IMPL_PROPOSED_NAMESPACE=experimental")
