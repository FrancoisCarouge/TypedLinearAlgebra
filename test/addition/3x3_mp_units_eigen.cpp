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
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

namespace {
//! @test Verifies the matrix addition operator of matrices whose rows and
//! columns repeat their element types, as the elements of a covariance repeat
//! over the points of a grid.
[[maybe_unused]] const auto test{[] -> int {
  using position = quantity<mp_units::isq::length[m]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;
  using dimensionless = quantity<mp_units::dimensionless[one]>;
  using rate = quantity<one / mp_units::isq::duration[s]>;
  using row_indexes = std::tuple<position, velocity, position>;
  using column_indexes = std::tuple<dimensionless, rate, dimensionless>;

  matrix<representation, row_indexes, column_indexes> lhs;
  matrix<representation, row_indexes, column_indexes> rhs;

  lhs.at<0, 0>(1. * m);
  lhs.at<0, 1>(2. * m / s);
  lhs.at<0, 2>(3. * m);
  lhs.at<1, 0>(4. * m / s);
  lhs.at<1, 1>(5. * m / s2);
  lhs.at<1, 2>(6. * m / s);
  lhs.at<2, 0>(7. * m);
  lhs.at<2, 1>(8. * m / s);
  lhs.at<2, 2>(9. * m);

  rhs.at<0, 0>(10. * m);
  rhs.at<0, 1>(20. * m / s);
  rhs.at<0, 2>(30. * m);
  rhs.at<1, 0>(40. * m / s);
  rhs.at<1, 1>(50. * m / s2);
  rhs.at<1, 2>(60. * m / s);
  rhs.at<2, 0>(70. * m);
  rhs.at<2, 1>(80. * m / s);
  rhs.at<2, 2>(90. * m);

  const matrix<representation, row_indexes, column_indexes> sum{lhs + rhs};

  assert((11. * m == sum.at<0, 0>()));
  assert((22. * m / s == sum.at<0, 1>()));
  assert((33. * m == sum.at<0, 2>()));
  assert((44. * m / s == sum.at<1, 0>()));
  assert((55. * m / s2 == sum.at<1, 1>()));
  assert((66. * m / s == sum.at<1, 2>()));
  assert((77. * m == sum.at<2, 0>()));
  assert((88. * m / s == sum.at<2, 1>()));
  assert((99. * m == sum.at<2, 2>()));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
