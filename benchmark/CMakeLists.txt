#[[ Typed Linear Algebra
Version 0.2.0
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

FetchContent_Declare(
  nanobench
  GIT_REPOSITORY "https://github.com/martinus/nanobench.git"
  GIT_SHALLOW TRUE
  FIND_PACKAGE_ARGS NAMES nanobench)
FetchContent_MakeAvailable(nanobench)

bench("eigen_product" "1" BACKENDS "eigen")

if(BUILD_BENCHMARKING)
  bench("eigen_product" "2" BACKENDS "eigen")
  bench("eigen_product" "4" BACKENDS "eigen")
  bench("eigen_product" "8" BACKENDS "eigen")
  bench("eigen_product" "16" BACKENDS "eigen")
  bench("eigen_product" "32" BACKENDS "eigen")
  bench("eigen_product" "64" BACKENDS "eigen")
  bench("eigen_product" "128" BACKENDS "eigen")
endif()

bench("typed_eigen_product" "1" BACKENDS "eigexed")

if(BUILD_BENCHMARKING)
  bench("typed_eigen_product" "2" BACKENDS "eigexed")
  bench("typed_eigen_product" "4" BACKENDS "eigexed")
  bench("typed_eigen_product" "8" BACKENDS "eigexed")
  bench("typed_eigen_product" "16" BACKENDS "eigexed")
  bench("typed_eigen_product" "32" BACKENDS "eigexed")
endif()

bench("mdspan_product" "1" BACKENDS "kokkos")

if(BUILD_BENCHMARKING)
  bench("mdspan_product" "2" BACKENDS "kokkos")
  bench("mdspan_product" "4" BACKENDS "kokkos")
  bench("mdspan_product" "8" BACKENDS "kokkos")
  bench("mdspan_product" "16" BACKENDS "kokkos")
  bench("mdspan_product" "32" BACKENDS "kokkos")
  bench("mdspan_product" "64" BACKENDS "kokkos")
  bench("mdspan_product" "128" BACKENDS "kokkos")
endif()

bench("mdspan_tuple_product" "1" BACKENDS "unit_std")

if(BUILD_BENCHMARKING)
  bench("mdspan_tuple_product" "2" BACKENDS "unit_std")
  bench("mdspan_tuple_product" "4" BACKENDS "unit_std")
  bench("mdspan_tuple_product" "8" BACKENDS "unit_std")
  bench("mdspan_tuple_product" "16" BACKENDS "unit_std")
endif()

bench("typed_mdspan_tupled_product" "1" BACKENDS "unit_std")

if(BUILD_BENCHMARKING)
  bench("typed_mdspan_tupled_product" "2" BACKENDS "unit_std")
  bench("typed_mdspan_tupled_product" "4" BACKENDS "unit_std")
  bench("typed_mdspan_tupled_product" "8" BACKENDS "unit_std")
  bench("typed_mdspan_tupled_product" "16" BACKENDS "unit_std")
endif()
