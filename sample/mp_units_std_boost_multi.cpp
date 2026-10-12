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

//! @file
//! @brief Unit safe linear algebra with mp-units over Boost.Multi memory.
//!
//! @details Demonstrate a variety of linear algebra operations with mp-units,
//! std::linalg, and std::mdspan over Boost.Multi memory. The application owns
//! the memory of every matrix in Boost.Multi two-dimension arrays. This library
//! composes a std::mdspan view of that memory as the matrix' linear algebra
//! backend with index typed as mp-units types. The arrays are
//! `boost::multi::inplace_array`, keeping their bounded memory in place, on the
//! stack: the sample never allocates and runs in the real-time context. A
//! nested braced list would allocate its rows: each array is constructed from
//! its extents, then assigned its elements in row-major order. This sample
//! explicitly uses double precision floating point numbers.

#include "fcarouge/linalg.hpp"

#include <boost/multi/array.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <mdspan>
#include <memory>
#include <tuple>

namespace fcarouge::sample {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

using position = quantity<mp_units::isq::length[m]>;
using velocity = quantity<mp_units::isq::velocity[m / s]>;
using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

template <typename RowIndexes, typename ColumnIndexes>
using matrix = matrix<representation, RowIndexes, ColumnIndexes>;

template <typename... Types>
using column_vector = column_vector<representation, Types...>;

template <typename... Types>
using row_vector = row_vector<representation, Types...>;

using state = column_vector<position, velocity, acceleration>;
using state_transpose = row_vector<position, velocity, acceleration>;
using estimate_uncertainty =
    matrix<std::tuple<position, velocity, acceleration>,
           std::tuple<position, velocity, acceleration>>;

using column_extents = std::extents<std::size_t, 3, 1>;
using row_extents = std::extents<std::size_t, 1, 3>;
using matrix_extents = std::extents<std::size_t, 3, 3>;

using storage = boost::multi::inplace_array<representation, 2>;

//! @brief Strongly typed linear algebra samples.
//!
//! @details A variety of activities of strongly typed linear algebra with
//! Boost.Multi, std::mdspan, std::linalg, and mp-units.
[[maybe_unused]] const auto sample{[] -> int {
  // The application owns the memory, here initialized in its representation.
  storage a0({3, 1}, 0.);
  a0.elements() = {3., 2., 1.};
  state x0{std::mdspan{std::to_address(a0.data_elements()), column_extents{}}};

  // Typed element access of the owned memory.
  assert(x0.at<0>() == 3. * m);
  assert(x0.at<1>() == 2. * m / s);
  assert(x0.at<2>() == 1. * m / s2);

  // Typed element assignment into the owned memory.
  x0.at<1>(2.5 * m / s);
  assert(std::ranges::equal(a0.elements(), std::array{3., 2.5, 1.}));

  // Multiplication with a scalar factor, in place.
  scale(3., x0);
  assert(std::ranges::equal(a0.elements(), std::array{9., 7.5, 3.}));

  storage a1({3, 1}, 0.);
  state x1{std::mdspan{std::to_address(a1.data_elements()), column_extents{}}};

  // Additions of two vectors of the same types.
  add(x0, x0, x1);
  assert(x1.at<0>() == 18. * m);
  assert(x1.at<1>() == 15. * m / s);
  assert(x1.at<2>() == 6. * m / s2);

  storage a2({1, 3}, 0.);
  a2.elements() = {3., 2., 1.};
  state_transpose xt2{
      std::mdspan{std::to_address(a2.data_elements()), row_extents{}}};

  storage a3({3, 3}, 0.);
  estimate_uncertainty p3{
      std::mdspan{std::to_address(a3.data_elements()), matrix_extents{}}};

  // Column-row matrix product, stored in row-major order.
  matrix_product(x0, xt2, p3);
  assert((p3.at<0, 0>() == 27. * m2));
  assert((p3.at<1, 0>() == 22.5 * m2 / s));
  assert((p3.at<2, 2>() == 3. * m2 / (s2 * s2)));
  assert(std::ranges::equal(
      a3.elements(), std::array{27., 18., 9., 22.5, 15., 7.5, 9., 6., 3.}));

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
