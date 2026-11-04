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
#include <tuple>
#include <units/length.h>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the conversions between sparse typed matrices of
//! nholthaus/units quantities of one dimension in different units, meters and
//! kilometers, keep the values of the stored elements and the structural zeros,
//! whatever the pattern.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using small_unit = units::length::meters<representation>;
  using large_unit = units::length::kilometers<representation>;
  using scalars = std::tuple<representation, representation>;

  constexpr auto top{sparse::make_pattern<2, 1>({{0, 0}})};

  using small_column =
      sparse_column_vector<representation, top, small_unit, small_unit>;
  using large_column =
      sparse_column_vector<representation, top, large_unit, large_unit>;
  using small_dense_column =
      sparse_column_vector<representation, sparse::dense<2, 1>, small_unit,
                           small_unit>;
  using small_matrix =
      sparse_matrix<representation, sparse::diagonal<2>,
                    std::tuple<small_unit, small_unit>, scalars>;
  using large_matrix =
      sparse_matrix<representation, sparse::diagonal<2>,
                    std::tuple<large_unit, large_unit>, scalars>;

  // A vector with a structural zero, on construction and on assignment.
  large_column large_operand;

  large_operand.at<0>(large_unit{1.});

  const small_column converted{large_operand};
  small_column assigned;

  assigned = large_operand;

  assert(converted.at<0>() == small_unit{1000.});
  assert(converted.at<1>() == small_unit{0.});
  assert(assigned.at<0>() == small_unit{1000.});
  assert(assigned.at<1>() == small_unit{0.});

  // A matrix with structural zeros.
  large_matrix large_diagonal;

  large_diagonal.at<0, 0>(large_unit{1.});
  large_diagonal.at<1, 1>(large_unit{2.});

  const small_matrix small_diagonal{large_diagonal};

  assert((small_diagonal.at<0, 0>() == small_unit{1000.}));
  assert((small_diagonal.at<0, 1>() == small_unit{0.}));
  assert((small_diagonal.at<1, 1>() == small_unit{2000.}));

  // Into a pattern storing more positions, never fewer.
  const small_dense_column widened{large_operand};

  assert(widened.at<0>() == small_unit{1000.});
  assert(widened.at<1>() == small_unit{0.});
  static_assert(not std::convertible_to<
                sparse_column_vector<representation, sparse::dense<2, 1>,
                                     large_unit, large_unit>,
                small_column>);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
