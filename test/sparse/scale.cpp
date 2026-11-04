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
//! @test Verifies the scaling, and the scalar division, of a matrix.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using tridiagonal_matrix = sparse_matrix<double, sparse::tridiagonal<3>>;

  const tridiagonal_matrix value{sparse::make_matrix<tridiagonal_matrix>(
      1., 2., sparse::structural_zero, 3., 4., 5., sparse::structural_zero, 6.,
      7.)};

  const tridiagonal_matrix right{value * 2.};
  const tridiagonal_matrix left{2. * value};
  const tridiagonal_matrix divided{value / 2.};

  assert((right.at<0, 0>() == 2.));
  assert((right.at<2, 2>() == 14.));
  assert((right.at<0, 2>() == 0.));
  assert((right == left));
  assert((divided.at<1, 2>() == 2.5));
  assert((divided.at<2, 0>() == 0.));
  assert((divided * 2. == value));

  // A factor of zero keeps the pattern.
  const tridiagonal_matrix none{value * 0.};

  assert((none == tridiagonal_matrix{}));

  // A singleton and a vector.
  const sparse_matrix<double, sparse::make_pattern<3, 1>({{1, 0}})> column{
      sparse::make_matrix<
          sparse_matrix<double, sparse::make_pattern<3, 1>({{1, 0}})>>(
          sparse::structural_zero, 4., sparse::structural_zero)};
  const auto scaled_column{column * 3.};

  assert(scaled_column.at<0>() == 0.);
  assert(scaled_column.at<1>() == 12.);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
