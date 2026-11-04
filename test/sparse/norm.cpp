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
//! @test Verifies the norms. The operands discriminate a wrong implementation:
//! negative elements for an absolute value, a non-square matrix for a row and
//! column mix-up, a row vector for a vector norm special case, a position
//! shared between rows and columns.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_pattern;

  // The largest absolute column sum is not the largest absolute row sum.
  constexpr auto wide{make_pattern<2, 3>({{0, 0}, {0, 1}, {1, 0}, {1, 2}})};

  const sparse_matrix<double, wide> value{
      sparse::make_matrix<sparse_matrix<double, wide>>(
          1., -2., sparse::structural_zero, -3., sparse::structural_zero, 4.)};

  assert(matrix_one_norm(value) == 4.);
  assert(matrix_inf_norm(value) == 7.);
  assert(matrix_frob_norm(value) == sparse::square_root(30.));

  // A row vector: the one norm is its largest absolute element, the infinity
  // norm is the sum of its absolute elements.
  const sparse_matrix<double, make_pattern<1, 3>({{0, 0}, {0, 2}})> row{
      sparse::make_matrix<
          sparse_matrix<double, make_pattern<1, 3>({{0, 0}, {0, 2}})>>(
          -3., sparse::structural_zero, 4.)};

  assert(matrix_one_norm(row) == 4.);
  assert(matrix_inf_norm(row) == 7.);
  assert(matrix_frob_norm(row) == 5.);
  assert(magnitude(row) == 5.);

  // A column vector, the other way around.
  const sparse_matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})> column{
      sparse::make_matrix<
          sparse_matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})>>(
          -3., sparse::structural_zero, 4.)};

  assert(matrix_one_norm(column) == 7.);
  assert(matrix_inf_norm(column) == 4.);
  assert(matrix_frob_norm(column) == 5.);
  assert(magnitude(column) == 5.);

  // A row or a column without any stored element.
  const sparse_matrix<double, make_pattern<2, 2>({{0, 0}})> corner{
      sparse::make_matrix<sparse_matrix<double, make_pattern<2, 2>({{0, 0}})>>(
          -2., sparse::structural_zero, sparse::structural_zero,
          sparse::structural_zero)};

  assert(matrix_one_norm(corner) == 2.);
  assert(matrix_inf_norm(corner) == 2.);
  assert(matrix_frob_norm(corner) == 2.);

  // Every negative column: the maximum is not started from a zero absent from
  // the absolute values, nor the first line.
  const sparse_matrix<double, make_pattern<2, 2>({{0, 1}, {1, 1}})> last_column{
      sparse::make_matrix<
          sparse_matrix<double, make_pattern<2, 2>({{0, 1}, {1, 1}})>>(
          sparse::structural_zero, -1., sparse::structural_zero, -2.)};

  assert(matrix_one_norm(last_column) == 3.);
  assert(matrix_inf_norm(last_column) == 2.);

  // The zero matrix.
  const sparse_matrix<double, sparse::zero<2, 2>> none;

  assert(matrix_one_norm(none) == 0.);
  assert(matrix_inf_norm(none) == 0.);
  assert(matrix_frob_norm(none) == 0.);

  // An integer one norm, needing no square root.
  const sparse_matrix<int, make_pattern<2, 2>({{0, 0}, {1, 0}})> integers{
      sparse::make_matrix<
          sparse_matrix<int, make_pattern<2, 2>({{0, 0}, {1, 0}})>>(
          -2, sparse::structural_zero, 3, sparse::structural_zero)};

  assert(matrix_one_norm(integers) == 5);
  assert(matrix_inf_norm(integers) == 3);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
