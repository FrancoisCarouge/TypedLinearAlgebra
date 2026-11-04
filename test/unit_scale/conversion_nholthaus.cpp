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
of the public at large_unit and to the detriment of our heirs and
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

#include <units/length.h>

#include <cassert>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the conversions between typed matrices of nholthaus/units
//! quantities of one dimension in different units, meters and kilometers, keep
//! the values, as the elements themselves do on conversion.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = units::length::meters<representation>;
  using large_unit = units::length::kilometers<representation>;
  using scalars = std::tuple<representation, representation>;
  using small_column = column_vector<representation, small_unit, small_unit>;
  using large_column = column_vector<representation, large_unit, large_unit>;
  using small_singleton = column_vector<representation, small_unit>;
  using large_singleton = column_vector<representation, large_unit>;
  using small_matrix =
      matrix<representation, std::tuple<small_unit, small_unit>, scalars>;
  using large_matrix =
      matrix<representation, std::tuple<large_unit, large_unit>, scalars>;

  const large_unit value{3.};
  const large_column source{value, value};

  // Copy construction keeps the value, whichever way the element converts.
  const small_column constructed{source};

  assert(constructed.at<0>() == value);
  assert(constructed.at<1>() == value);
  assert(constructed == source);

  // Copy initialization, the implicit conversion.
  const small_column initialized = source;

  assert(initialized.at<0>() == value);

  // Move construction.
  const small_column moved{large_column{value, value}};

  assert(moved.at<0>() == value);

  // Copy assignment.
  small_column assigned{value, value};

  assigned = source;

  assert(assigned.at<0>() == value);

  // Move assignment.
  small_column move_assigned{value, value};

  move_assigned = large_column{value, value};

  assert(move_assigned.at<0>() == value);

  // From an expression template.
  const small_column from_expression{source + source};

  assert(from_expression.at<0>() == value + value);

  assigned = source + source;

  assert(assigned.at<1>() == value + value);

  // A singleton, whose addition already honors the units.
  const large_singleton singleton_source{value};
  const small_singleton singleton{singleton_source};

  assert(singleton.at<>() == value);

  // A matrix.
  large_matrix matrix_source;

  matrix_source.at<0, 1>(value);

  const small_matrix matrix_converted{matrix_source};

  assert((matrix_converted.at<0, 1>() == value));
  assert(matrix_converted == matrix_source);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
