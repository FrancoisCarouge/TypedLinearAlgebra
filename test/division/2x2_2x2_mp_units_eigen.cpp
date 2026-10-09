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
//! @test Verifies the matrix by matrix division operator when the operands'
//! first column indexes carry different units: dimensionless columns for the
//! dividend, length columns for the divisor. The shape of a Kalman filter gain
//! `P * transposed(H) / S`, where the indexes of the quotient must not be
//! derived as if both first columns were alike.
[[maybe_unused]] const auto test{[] -> int {
  using dimensionless = quantity<one>;
  using length = quantity<mp_units::isq::length[m]>;
  using area = quantity<mp_units::isq::area[m2]>;
  using area_rate =
      quantity<mp_units::isq::area[m2] / mp_units::isq::duration[s]>;
  using lhs_row_indexes = std::tuple<area, area_rate>;
  using lhs_column_indexes = std::tuple<dimensionless, dimensionless>;
  using rhs_indexes = std::tuple<length, length>;

  matrix<representation, lhs_row_indexes, lhs_column_indexes> lhs;
  matrix<representation, rhs_indexes, rhs_indexes> rhs;

  lhs.at<0, 0>(6. * m2);
  lhs.at<0, 1>(5. * m2);
  lhs.at<1, 0>(11. * m2 / s);
  lhs.at<1, 1>(1. * m2 / s);

  rhs.at<0, 0>(4. * m2);
  rhs.at<0, 1>(1. * m2);
  rhs.at<1, 0>(1. * m2);
  rhs.at<1, 1>(2. * m2);

  const auto quotient{lhs / rhs};
  const auto product{quotient * rhs};

  assert(abs(quotient.at<0, 0>() - 1. * one) < 1e-9 * one);
  assert(abs(quotient.at<0, 1>() - 2. * one) < 1e-9 * one);
  assert(abs(quotient.at<1, 0>() - 3. / s) < 1e-9 / s);
  assert(abs(quotient.at<1, 1>() + 1. / s) < 1e-9 / s);
  assert(abs(product.at<0, 0>() - 6. * m2) < 1e-9 * m2);
  assert(abs(product.at<0, 1>() - 5. * m2) < 1e-9 * m2);
  assert(abs(product.at<1, 0>() - 11. * m2 / s) < 1e-9 * m2 / s);
  assert(abs(product.at<1, 1>() - 1. * m2 / s) < 1e-9 * m2 / s);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
