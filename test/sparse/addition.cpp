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
#include <concepts>
#include <type_traits>

namespace fcarouge::test {
namespace {
//! @test Verifies the addition of matrices of overlapping and disjoint
//! patterns.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_pattern;

  constexpr auto first_row{make_pattern<2, 2>({{0, 0}, {0, 1}})};
  constexpr auto first_column{make_pattern<2, 2>({{0, 0}, {1, 0}})};
  constexpr auto opposite{make_pattern<2, 2>({{1, 1}})};

  const sparse_matrix<double, first_row> lhs{
      sparse::make_matrix<sparse_matrix<double, first_row>>(
          1., 2., sparse::structural_zero, sparse::structural_zero)};
  const sparse_matrix<double, first_column> rhs{
      sparse::make_matrix<sparse_matrix<double, first_column>>(
          3., sparse::structural_zero, 4., sparse::structural_zero)};
  const sparse_matrix<double, opposite> other{
      sparse::make_matrix<sparse_matrix<double, opposite>>(
          sparse::structural_zero, sparse::structural_zero,
          sparse::structural_zero, 5.)};

  // Overlapping: the union of the patterns.
  const auto sum{lhs + rhs};

  static_assert(
      std::same_as<std::remove_cvref_t<decltype(sum)>,
                   sparse_matrix<double, sparse::pattern_union<first_row,
                                                               first_column>>>);
  assert((sum.at<0, 0>() == 4.));
  assert((sum.at<0, 1>() == 2.));
  assert((sum.at<1, 0>() == 4.));
  assert((sum.at<1, 1>() == 0.));

  // Disjoint.
  const auto disjoint{lhs + other};

  assert((disjoint.at<0, 0>() == 1.));
  assert((disjoint.at<0, 1>() == 2.));
  assert((disjoint.at<1, 0>() == 0.));
  assert((disjoint.at<1, 1>() == 5.));

  // Commutative.
  assert((lhs + rhs == rhs + lhs));

  // A zero matrix is neutral, and the pattern does not grow.
  const auto neutral{lhs + sparse_matrix<double, sparse::zero<2, 2>>{}};

  static_assert(std::same_as<std::remove_cvref_t<decltype(neutral)>,
                             sparse_matrix<double, first_row>>);
  assert((neutral == lhs));

  // Elements cancelling out at run time stay stored: a pattern is structural.
  const auto cancelled{lhs + (-lhs)};

  static_assert(std::same_as<std::remove_cvref_t<decltype(cancelled)>,
                             sparse_matrix<double, first_row>>);
  assert((cancelled.at<0, 0>() == 0.));

  // Run-time values reach the result, through a dense reference.
  const sparse_matrix<double, sparse::dense<2, 2>> reference{{4., 2.},
                                                             {4., 0.}};

  assert((sum == reference));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
