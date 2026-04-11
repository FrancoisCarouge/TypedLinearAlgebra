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

add_executable(typed_linear_algebra_benchmark_matrix_product_mdspan_driver
               "matrix_product_mdspan.cpp")
target_link_libraries(
  typed_linear_algebra_benchmark_matrix_product_mdspan_driver
  PRIVATE typed_linear_algebra_options typed_linear_algebra_kokkos
          nanobench::nanobench)
add_test(NAME typed_linear_algebra_benchmark_matrix_product_mdspan
         COMMAND typed_linear_algebra_benchmark_matrix_product_mdspan_driver)

add_executable(typed_linear_algebra_benchmark_matrix_product_tuple_driver
               "matrix_product_tuple.cpp")
target_link_libraries(
  typed_linear_algebra_benchmark_matrix_product_tuple_driver
  PRIVATE typed_linear_algebra_options typed_linear_algebra_kokkos
          nanobench::nanobench)
add_test(NAME typed_linear_algebra_benchmark_matrix_product_tuple
         COMMAND typed_linear_algebra_benchmark_matrix_product_tuple_driver)

         add_executable(typed_linear_algebra_benchmark_eigen_product_driver
               "eigen_product.cpp")
target_link_libraries(
  typed_linear_algebra_benchmark_eigen_product_driver
  PRIVATE typed_linear_algebra_options typed_linear_algebra_eigen
          nanobench::nanobench)
add_test(NAME typed_linear_algebra_benchmark_eigen_product
         COMMAND typed_linear_algebra_benchmark_eigen_product_driver)

#   add_executable(typed_linear_algebra_benchmark_eigen_product_typed_driver
# "eigen_product_typed.cpp")
# target_link_libraries(
#   typed_linear_algebra_benchmark_eigen_product_typed_driver
#   PRIVATE typed_linear_algebra_options typed_linear_algebra_eigexed
#           nanobench::nanobench)
# add_test(NAME typed_linear_algebra_benchmark_eigen_product_typed
#          COMMAND typed_linear_algebra_benchmark_eigen_product_typed_driver)
