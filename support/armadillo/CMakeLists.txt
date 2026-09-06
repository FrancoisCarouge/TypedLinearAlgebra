#[[ Typed Linear Algebra
Version 0.3.0
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
  armadillo
  GIT_REPOSITORY "https://gitlab.com/conradsnicta/armadillo-code.git"
  GIT_TAG "15.6.0"
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR "include" SYSTEM FIND_PACKAGE_ARGS NAMES Armadillo)
FetchContent_MakeAvailable(armadillo)

add_library(typed_linear_algebra_armadillo INTERFACE)
target_sources(
  typed_linear_algebra_armadillo
  INTERFACE FILE_SET
            "typed_linear_algebra_headers"
            TYPE
            "HEADERS"
            FILES
            "fcarouge/armadillo.hpp"
            "fcarouge/linalg.hpp")
target_link_libraries(typed_linear_algebra_armadillo INTERFACE tlinalg)
target_compile_definitions(
  typed_linear_algebra_armadillo
  INTERFACE "ARMA_DONT_USE_WRAPPER"
            "ARMA_DONT_USE_LAPACK"
            "ARMA_DONT_USE_BLAS"
            "ARMA_DONT_USE_ARPACK"
            "ARMA_DONT_USE_SUPERLU"
            "ARMA_DONT_USE_HDF5"
            "ARMA_DONT_USE_FFTW3")

# if(TARGET armadillo)
target_link_libraries(typed_linear_algebra_armadillo INTERFACE armadillo)
# else() target_include_directories(typed_linear_algebra_armadillo SYSTEM
# INTERFACE "${armadillo_SOURCE_DIR}/include") endif()
