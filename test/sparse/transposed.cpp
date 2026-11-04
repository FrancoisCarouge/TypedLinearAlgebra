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
//! @test Verifies the transposition moves each element to its mirrored
//! position, on non-square and unsymmetric patterns.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_pattern;

  constexpr auto wide{make_pattern<2, 3>({{0, 0}, {0, 2}, {1, 1}})};

  const sparse_matrix<double, wide> value{
      sparse::make_matrix<sparse_matrix<double, wide>>(
          1., sparse::structural_zero, 2., sparse::structural_zero, 3.,
          sparse::structural_zero)};
  const auto mirrored{transposed(value)};

  static_assert(decltype(mirrored)::rows == 3 and
                decltype(mirrored)::columns == 2);
  assert((mirrored.at<0, 0>() == 1.));
  assert((mirrored.at<2, 0>() == 2.));
  assert((mirrored.at<1, 1>() == 3.));
  assert((mirrored.at<0, 1>() == 0.));
  assert((mirrored.at<1, 0>() == 0.));
  assert((mirrored.at<2, 1>() == 0.));

  assert((transposed(mirrored) == value));

  // The stored order changes: a column then row reading of the storage.
  const sparse_matrix<double, sparse::lower_triangular<3>> lower{
      sparse::make_matrix<sparse_matrix<double, sparse::lower_triangular<3>>>(
          1., sparse::structural_zero, sparse::structural_zero, 2., 3.,
          sparse::structural_zero, 4., 5., 6.)};
  const auto upper{transposed(lower)};

  assert((upper.at<0, 0>() == 1.));
  assert((upper.at<0, 1>() == 2.));
  assert((upper.at<0, 2>() == 4.));
  assert((upper.at<1, 1>() == 3.));
  assert((upper.at<1, 2>() == 5.));
  assert((upper.at<2, 2>() == 6.));
  assert((upper.at<1, 0>() == 0.));

  // A row vector to a column vector.
  const sparse_matrix<double, make_pattern<1, 3>({{0, 1}})> row{
      sparse::make_matrix<sparse_matrix<double, make_pattern<1, 3>({{0, 1}})>>(
          sparse::structural_zero, 7., sparse::structural_zero)};
  const auto column{transposed(row)};

  assert(column.at<0>() == 0.);
  assert(column.at<1>() == 7.);
  assert(column.at<2>() == 0.);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
