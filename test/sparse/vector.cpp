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
//! @test Verifies a sparse vector with heterogeneous element types, its
//! constructions, and its decomposition.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_pattern;

  constexpr auto ends{make_pattern<3, 1>({{0, 0}, {2, 0}})};

  using value_type = sparse_column_vector<double, ends, double, double, double>;

  // A structural zero is explicit in the construction of the elements.
  const value_type value{
      sparse::make_matrix<value_type>(1., sparse::structural_zero, 3.)};

  assert(value.at<0>() == 1.);
  assert(value.at<1>() == 0.);
  assert(value.at<2>() == 3.);

  // A structured binding reads the structural zero by value.
  const auto [first, second, third]{value};

  assert(first == 1.);
  assert(second == 0.);
  assert(third == 3.);

  // A row vector.
  using row_type = sparse_row_vector<double, make_pattern<1, 3>({{0, 1}}),
                                     double, double, double>;

  const row_type row{sparse::make_matrix<row_type>(sparse::structural_zero, 2.,
                                                   sparse::structural_zero)};

  assert(row.at<0>() == 0.);
  assert(row.at<1>() == 2.);

  const auto [row_first, row_second, row_third]{row};

  assert(row_first == 0.);
  assert(row_second == 2.);
  assert(row_third == 0.);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
