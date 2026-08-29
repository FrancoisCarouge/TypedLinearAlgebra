/* Typed Linear Algebra
Version 0.3.0
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

//! @file
//! @brief Unit safe linear algebra with Au and Eigen.
//!
//! @details Demonstrate a variety of linear algebra operations with Au and
//! Eigen. This library composes Eigen as the matrix' linear algebra backend
//! with indexes typed as Au types. This sample explicitly uses double
//! precision floating point numbers. This sample uses Eigen linear algebra as
//! the linear algebra backend. This sample uses Au types for the strongly
//! typed units.

#include "fcarouge/linalg.hpp"

#include <au/std_format.hh>
#include <au/units/meters.hh>
#include <au/units/seconds.hh>

#include <cassert>
#include <format>
#include <print>
#include <tuple>

namespace fcarouge::sample {
namespace {
using representation = double;
using position = au::QuantityD<au::Meters>;
using velocity = au::QuantityD<au::UnitQuotientT<au::Meters, au::Seconds>>;
using acceleration = au::QuantityD<
    au::UnitQuotientT<au::Meters, au::UnitPowerT<au::Seconds, 2>>>;

template <typename RowIndexes, typename ColumnIndexes>
using matrix = matrix<representation, RowIndexes, ColumnIndexes>;

template <typename... Types>
using column_vector = column_vector<representation, Types...>;

template <typename... Types>
using row_vector = row_vector<representation, Types...>;

//! @brief Strongly typed linear algebra samples.
//!
//! @details A variety of activities of strongly typed linear algebra with
//! Eigen and Au.
[[maybe_unused]] const auto sample{[] -> int {
  using au::symbols::m;
  using au::symbols::s;

  constexpr auto m2{au::squared(m)};
  constexpr auto s2{au::squared(s)};

  // Set up a heterogenous column vector type for the sample.
  using state = column_vector<position, velocity, acceleration>;

  // A vector of quantities:
  state x0{3. * m, 2. * m / s, 1. * m / s2};

  // Printable.
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[3 m], [2 m / s], [1 m / s^2]]");

  // Element assignment and access.
  x0.at<1>(2.5 * m / s);
  auto x0_1{x0.at<1>()};
  assert(x0_1 == 2.5 * m / s);
  assert(std::format("{}", x0_1) == "2.5 m / s");

  // Multiplication with a scalar factor.
  state x1{x0 * 3.};
  assert(std::format("{}", x1) == "[[9 m], [7.5 m / s], [3 m / s^2]]");

  // Division with a scalar divisor.
  state x2{x1 / 2.};
  assert(std::format("{}", x2) == "[[4.5 m], [3.75 m / s], [1.5 m / s^2]]");

  // Substraction of two vectors of the same types.
  state x3{x2 - x0};
  assert(std::format("{}", x3) == "[[1.5 m], [1.25 m / s], [0.5 m / s^2]]");

  // Additions of two vectors of the same types.
  state x4{x3 + x3};
  assert(std::format("{}", x4) == "[[3 m], [2.5 m / s], [1 m / s^2]]");

  // Matrix and uniform typed access.
  using position_2d_uncertainty =
      matrix<std::tuple<position, position>, std::tuple<position, position>>;
  position_2d_uncertainty p0;

  p0.at<0, 1>(9. * m2);
  assert((p0[0, 1] == 9. * m2));

  p0.at<0, 1>(16. * m2);
  assert((p0(0, 1) == 16. * m2));

  // Matrix multiplication: column vector by row vector outer product,
  // producing an area typed matrix.
  const column_vector<position, position> c{1. * m, 2. * m};
  const row_vector<position, position> r{3. * m, 4. * m};
  const matrix<std::tuple<position, position>, std::tuple<position, position>>
      cr{c * r};

  assert((cr.at<0, 0>() == 3. * m2));
  assert((cr.at<0, 1>() == 4. * m2));
  assert((cr.at<1, 0>() == 6. * m2));
  assert((cr.at<1, 1>() == 8. * m2));

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
