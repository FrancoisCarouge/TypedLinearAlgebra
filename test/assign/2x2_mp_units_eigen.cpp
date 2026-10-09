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
#include <tuple>

namespace fcarouge::test {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::one;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;
using mp_units::si::unit_symbols::s;

namespace {
//! @test Verifies the assignment operator between matrices whose element types
//! agree while their indexes split them differently: all units in the rows and
//! dimensionless columns for the source, units in both for the destination.
[[maybe_unused]] const auto test{[] -> int {
  using dimensionless = quantity<one>;
  using length = quantity<mp_units::isq::length[m]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;
  using area = quantity<mp_units::isq::length[m] * mp_units::isq::length[m]>;
  using area_rate =
      quantity<mp_units::isq::velocity[m / s] * mp_units::isq::length[m]>;
  using source_row_indexes = std::tuple<area, area_rate>;
  using source_column_indexes = std::tuple<dimensionless, dimensionless>;
  using destination_row_indexes = std::tuple<length, velocity>;
  using destination_column_indexes = std::tuple<length, length>;

  matrix<representation, source_row_indexes, source_column_indexes> source;
  matrix<representation, destination_row_indexes, destination_column_indexes>
      destination;

  source.at<0, 0>(1. * m2);
  source.at<0, 1>(2. * m2);
  source.at<1, 0>(3. * m2 / s);
  source.at<1, 1>(4. * m2 / s);

  destination = source;

  assert((destination.at<0, 0>() == 1. * m2));
  assert((destination.at<0, 1>() == 2. * m2));
  assert((destination.at<1, 0>() == 3. * m2 / s));
  assert((destination.at<1, 1>() == 4. * m2 / s));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
