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

# FetchContent_Declare( TBB GIT_REPOSITORY
# "https://github.com/uxlfoundation/oneTBB" GIT_SHALLOW TRUE FIND_PACKAGE_ARGS
# NAMES TBB) FetchContent_MakeAvailable(TBB)

# FetchContent_Declare( Kokkos GIT_REPOSITORY "https://github.com/kokkos/kokkos"
# GIT_TAG "develop" FIND_PACKAGE_ARGS NAMES Kokkos)
# FetchContent_MakeAvailable(Kokkos)

# FetchContent_Declare( KokkosKernels GIT_REPOSITORY
# "https://github.com/kokkos/kokkos-kernels" GIT_TAG "develop" FIND_PACKAGE_ARGS
# NAMES KokkosKernels) FetchContent_MakeAvailable(KokkosKernels)

# FetchContent_Declare( LinAlg GIT_REPOSITORY
# "https://github.com/kokkos/stdBLAS" GIT_TAG "main" FIND_PACKAGE_ARGS NAMES
# LinAlg) FetchContent_MakeAvailable(LinAlg)

FetchContent_Declare(
  mdspan
  GIT_REPOSITORY "https://github.com/kokkos/mdspan"
  GIT_TAG "stable"
  FIND_PACKAGE_ARGS NAMES mdspan)
FetchContent_MakeAvailable(mdspan)

add_library(typed_linear_algebra_kokkos INTERFACE)
# target_sources( typed_linear_algebra_kokkos INTERFACE FILE_SET
# "typed_linear_algebra_headers" TYPE "HEADERS" FILES "fcarouge/kokkos.hpp"
# "fcarouge/linalg.hpp")
target_link_libraries(
  typed_linear_algebra_kokkos INTERFACE # TBB::tbb Kokkos::kokkos LinAlg::linalg
                                        mdspan::mdspan tlinalg)
