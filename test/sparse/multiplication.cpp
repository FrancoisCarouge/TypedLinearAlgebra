/* Typed Linear Algebra
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

For more information, please refer to <https://unlicense.org> */

#include "fcarouge/linalg.hpp"

#include <cassert>

namespace fcarouge::test {
namespace {
//! @test Verifies the matrix multiplication across shapes and patterns, at run
//! time.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::diagonal;
  using sparse::make_pattern;
  using sparse::tridiagonal;

  const sparse_matrix<double, diagonal<3>> scale_matrix{
      sparse::make_matrix<sparse_matrix<double, diagonal<3>>>(
          2., sparse::structural_zero, sparse::structural_zero,
          sparse::structural_zero, 3., sparse::structural_zero,
          sparse::structural_zero, sparse::structural_zero, 4.)};
  const sparse_matrix<double, tridiagonal<3>> band{
      sparse::make_matrix<sparse_matrix<double, tridiagonal<3>>>(
          1., 2., sparse::structural_zero, 3., 4., 5., sparse::structural_zero,
          6., 7.)};

  // Left scales the rows, right scales the columns: not commutative.
  const auto left{scale_matrix * band};
  const auto right{band * scale_matrix};

  assert((left.at<0, 1>() == 4.));
  assert((left.at<1, 0>() == 9.));
  assert((left.at<2, 1>() == 24.));
  assert((right.at<0, 1>() == 6.));
  assert((right.at<1, 0>() == 6.));
  assert((right.at<2, 1>() == 18.));
  assert((not(left == right)));

  // A band squared reaches further: every term of the middle row.
  const auto square{band * band};

  assert((square.at<0, 0>() == 7.));
  assert((square.at<0, 1>() == 10.));
  assert((square.at<0, 2>() == 10.));
  assert((square.at<1, 1>() == 52.));
  assert((square.at<2, 0>() == 18.));
  assert((square.at<2, 2>() == 79.));

  // A non-square product.
  constexpr auto wide{make_pattern<2, 3>({{0, 0}, {0, 2}, {1, 1}})};
  constexpr auto tall{make_pattern<3, 2>({{0, 1}, {1, 0}, {2, 1}})};

  const sparse_matrix<double, wide> lhs{
      sparse::make_matrix<sparse_matrix<double, wide>>(
          1., sparse::structural_zero, 2., sparse::structural_zero, 3.,
          sparse::structural_zero)};
  const sparse_matrix<double, tall> rhs{
      sparse::make_matrix<sparse_matrix<double, tall>>(
          sparse::structural_zero, 4., 5., sparse::structural_zero,
          sparse::structural_zero, 6.)};

  const auto product{lhs * rhs};

  static_assert(decltype(product)::rows == 2 and
                decltype(product)::columns == 2);
  assert((product.at<0, 0>() == 0.));
  assert((product.at<0, 1>() == 16.));
  assert((product.at<1, 0>() == 15.));
  assert((product.at<1, 1>() == 0.));

  const auto reversed{rhs * lhs};

  static_assert(decltype(reversed)::rows == 3 and
                decltype(reversed)::columns == 3);
  assert((reversed.at<0, 1>() == 12.));
  assert((reversed.at<1, 0>() == 5.));
  assert((reversed.at<1, 2>() == 10.));
  assert((reversed.at<2, 1>() == 18.));
  assert((reversed.at<0, 0>() == 0.));

  // A matrix by a vector, a sparse vector stays sparse.
  constexpr auto first_element{make_pattern<3, 1>({{0, 0}})};

  const sparse_matrix<double, first_element> column{
      sparse::make_matrix<sparse_matrix<double, first_element>>(
          2., sparse::structural_zero, sparse::structural_zero)};
  const auto image{band * column};

  assert(image.at<0>() == 2.);
  assert(image.at<1>() == 6.);
  assert(image.at<2>() == 0.);
  static_assert(decltype(image)::matrix::structure.count == 2);

  // A row by a column is a singleton, a column by a row is a matrix.
  constexpr auto row_pattern{make_pattern<1, 3>({{0, 0}, {0, 2}})};
  constexpr auto column_pattern{make_pattern<3, 1>({{0, 0}, {1, 0}})};

  const sparse_matrix<double, row_pattern> row_vector_value{
      sparse::make_matrix<sparse_matrix<double, row_pattern>>(
          3., sparse::structural_zero, 4.)};
  const sparse_matrix<double, column_pattern> column_vector_value{
      sparse::make_matrix<sparse_matrix<double, column_pattern>>(
          5., 6., sparse::structural_zero)};

  [[maybe_unused]] const double inner{row_vector_value * column_vector_value};
  const auto outer{column_vector_value * row_vector_value};

  assert(inner == 15.);
  assert((outer.at<0, 0>() == 15.));
  assert((outer.at<0, 2>() == 20.));
  assert((outer.at<1, 2>() == 24.));
  assert((outer.at<2, 2>() == 0.));
  static_assert(decltype(outer)::matrix::structure.count == 4);

  // Associative.
  assert(((scale_matrix * band) * band == scale_matrix * (band * band)));
  assert((transposed(scale_matrix * band) ==
          transposed(band) * transposed(scale_matrix)));

  // Distributive, over a wider pattern than either term.
  assert((band * (scale_matrix + band) == band * scale_matrix + band * band));

  // A zero matrix absorbs.
  const auto none{band * sparse_matrix<double, sparse::zero<3, 3>>{}};

  static_assert(decltype(none)::matrix::structure.count == 0);
  assert((none.at<1, 1>() == 0.));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
