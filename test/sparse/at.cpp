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
#include <cstddef>

namespace fcarouge::test {
using literals::operator""_i;

namespace {
//! @test Verifies the element access of a stored element and a structural zero,
//! through each accessor form, at run time.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using diagonal_matrix = sparse_matrix<double, sparse::diagonal<3>>;

  diagonal_matrix value;

  assert((value.at<0, 0>() == 0.));
  assert((value.at<0, 1>() == 0.));

  value.at<0, 0>(1.);
  value.at<1, 1>(2.);
  value.at<2, 2>(3.);

  assert((value.at<0, 0>() == 1.));
  assert((value.at<1, 1>() == 2.));
  assert((value.at<2, 2>() == 3.));
  assert((value.at<0, 1>() == 0.));
  assert((value.at<2, 0>() == 0.));

  assert((value[1, 1] == 2.));
  assert((value[1, 2] == 0.));
  assert((value(2, 2) == 3.));
  assert((value(2, 1) == 0.));
  assert((value[1_i, 1_i] == 2.));
  assert((value[1_i, 2_i] == 0.));
  assert((value(2_i, 2_i) == 3.));

  // Run-time indexes reach every position, for a uniform matrix.
  double trace{0.};
  double others{0.};

  for (std::size_t row{0}; row < 3; ++row) {
    for (std::size_t column{0}; column < 3; ++column) {
      (row == column ? trace : others) += value[row, column];
    }
  }

  assert(trace == 6.);
  assert(others == 0.);

  // A read is by value, never a reference, const or not.
  const diagonal_matrix constant{value};

  static_assert(std::same_as<decltype(value.at<0, 0>()), double>);
  static_assert(std::same_as<decltype(value.at<0, 1>()), double>);
  static_assert(std::same_as<decltype(constant.at<0, 0>()), double>);
  static_assert(std::same_as<decltype(value[0, 0]), double>);
  static_assert(std::same_as<decltype(value(0, 0)), double>);

  // A vector's linear index.
  sparse_matrix<double, sparse::make_pattern<4, 1>({{1, 0}, {3, 0}})> column;

  column.at<1>(5.);
  column.at<3>(6.);

  assert(column.at<0>() == 0.);
  assert(column.at<1>() == 5.);
  assert(column.at<2>() == 0.);
  assert(column.at<3>() == 6.);

  sparse_matrix<double, sparse::make_pattern<1, 4>({{0, 2}})> row;

  row.at<2>(7.);

  assert(row.at<1>() == 0.);
  assert(row.at<2>() == 7.);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
