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

#include <functional>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the matrix product is only available for terms of one scale,
//! with meters and kilometers.
[[maybe_unused]] const auto test{[] -> int {
  using mp_units::si::unit_symbols::km;
  using mp_units::si::unit_symbols::m;

  using small_unit =
      mp_units::quantity<mp_units::isq::length[m], representation>;
  using large_unit =
      mp_units::quantity<mp_units::isq::length[km], representation>;
  using gain =
      matrix<representation, std::tuple<representation, representation>,
             std::tuple<std::identity, std::identity>>;
  using small_column = column_vector<representation, small_unit, small_unit>;
  using large_column = column_vector<representation, large_unit, large_unit>;
  using mixed_column = column_vector<representation, small_unit, large_unit>;
  using mixed_row = row_vector<representation, small_unit, large_unit>;
  using scalar_row = row_vector<representation, representation, representation>;
  using scalar_column =
      column_vector<representation, representation, representation>;
  using mixed_columns =
      matrix<representation, std::tuple<representation, representation>,
             std::tuple<small_unit, large_unit>>;

  // Terms of one scale.
  static_assert(multipliable<gain, small_column>);
  static_assert(multipliable<gain, large_column>);
  static_assert(multipliable<scalar_row, small_column>);

  // Terms of different scales in a sum, whichever operand holds the units.
  static_assert(not multipliable<gain, mixed_column>);
  static_assert(not multipliable<scalar_row, mixed_column>);
  static_assert(not multipliable<mixed_row, gain>);
  static_assert(not multipliable<mixed_columns, scalar_column>);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
