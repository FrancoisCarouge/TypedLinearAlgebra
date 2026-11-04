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
//! @test Verifies the dot product only reaches the elements stored by both
//! vectors, of any orientation.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_pattern;

  constexpr auto ends{make_pattern<4, 1>({{0, 0}, {3, 0}})};
  constexpr auto inner{make_pattern<4, 1>({{1, 0}, {2, 0}})};
  constexpr auto overlap{make_pattern<4, 1>({{0, 0}, {1, 0}, {3, 0}})};
  constexpr auto overlap_row{make_pattern<1, 4>({{0, 1}, {0, 3}})};

  const sparse_matrix<double, ends> lhs{
      sparse::make_matrix<sparse_matrix<double, ends>>(
          2., sparse::structural_zero, sparse::structural_zero, -3.)};
  const sparse_matrix<double, inner> disjoint{
      sparse::make_matrix<sparse_matrix<double, inner>>(
          sparse::structural_zero, 5., 7., sparse::structural_zero)};
  const sparse_matrix<double, overlap> rhs{
      sparse::make_matrix<sparse_matrix<double, overlap>>(
          4., 5., sparse::structural_zero, 6.)};
  const sparse_matrix<double, overlap_row> row{
      sparse::make_matrix<sparse_matrix<double, overlap_row>>(
          sparse::structural_zero, 5., sparse::structural_zero, 6.)};

  assert(dot(lhs, rhs) == -10.);
  assert(dot(rhs, lhs) == -10.);
  assert(dot(lhs, disjoint) == 0.);
  assert(dot(lhs, lhs) == 13.);
  assert(dot(rhs, row) == 61.);
  assert(dot(row, rhs) == 61.);
  assert(dot(lhs, row) == -18.);

  // Against the empty vector.
  const sparse_matrix<double, sparse::zero<4, 1>> none;

  assert(dot(lhs, none) == 0.);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
