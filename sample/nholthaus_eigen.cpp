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
//! @brief Unit safe linear algebra with nholthaus/units and Eigen.
//!
//! @details Demonstrate a variety of linear algebra operations with the
//! nholthaus/units library and Eigen. This library composes Eigen as the
//! matrix' linear algebra backend with indexes typed as `units::unit`
//! containers. This sample explicitly uses double precision floating point
//! numbers as the underlying representation.

#include "fcarouge/linalg.hpp"

#include <units/acceleration.h>
#include <units/area.h>
#include <units/length.h>
#include <units/velocity.h>

#include <cassert>
#include <tuple>

namespace fcarouge::sample {
namespace {
using units::m;
using units::m2;
using units::mps;
using units::mps2;

using representation = double;
using position = units::length::meters<representation>;
using velocity = units::velocity::meters_per_second<representation>;
using acceleration =
    units::acceleration::meters_per_second_squared<representation>;

template <typename RowIndexes, typename ColumnIndexes>
using matrix = matrix<representation, RowIndexes, ColumnIndexes>;

template <typename... Types>
using column_vector = column_vector<representation, Types...>;

template <typename... Types>
using row_vector = row_vector<representation, Types...>;

//! @brief Strongly typed linear algebra samples.
//!
//! @details A variety of activities of strongly typed linear algebra with
//! Eigen and nholthaus/units.
[[maybe_unused]] const auto sample{[] -> int {
  // Set up a heterogenous column vector type for the sample.
  using state = column_vector<position, velocity, acceleration>;

  // A vector of quantities.
  state x0{3. * m, 2. * mps, 1. * mps2};

  // Element assignment and access.
  x0.at<1>(2.5 * mps);
  assert(x0.at<1>() == 2.5 * mps);

  // Multiplication with a scalar factor.
  const state x1{x0 * 3.};
  assert(x1.at<0>() == 9. * m);
  assert(x1.at<1>() == 7.5 * mps);
  assert(x1.at<2>() == 3. * mps2);

  // Division with a scalar divisor.
  const state x2{x1 / 2.};
  assert(x2.at<0>() == 4.5 * m);

  // Substraction then addition of vectors of the same types.
  const state x3{x2 - x0};
  const state x4{x3 + x3};
  assert(x4.at<0>() == 3. * m);

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
