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
//! @brief Unit safe linear algebra with mp-units and Eigen.
//!
//! @details Demonstrate a variety of linear algebra operations with mp-units
//! and Eigen. This library composes Eigen as the matrix' linear algebra backend
//! with index typed as mp-units types. This sample explicitly uses double
//! precision floating point numbers. This sample uses Eigen linear algebra as
//! the linear algebra backend. This sample uses mp-units types for the strongly
//! typed units.

#include "fcarouge/linalg.hpp"

#include <cstddef>
#include <format>
#include <linalg>
#include <mdspan>
#include <print>
#include <tuple>
#include <type_traits>
#include <vector>

namespace fcarouge::sample {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

template <typename RowIndexes, typename ColumnIndexes>
using matrix = matrix<representation, RowIndexes, ColumnIndexes>;

template <typename... Types>
using column_vector = column_vector<representation, Types...>;

template <typename... Types>
using row_vector = row_vector<representation, Types...>;

template <std::size_t Rows>
using column_extents = std::extents<std::size_t, Rows, 1>;

template <std::size_t Columns>
using row_extents = std::extents<std::size_t, 1, Columns>;

template <typename Extents>
constexpr std::size_t extents_size{[] {
  std::size_t size{1};
  Extents extents{};
  for (std::size_t i{0}; i < Extents::rank(); ++i) {
    size *= extents.extent(i);
  }
  return size;
}()};

//! @brief Strongly typed linear algebra samples.
//!
//! @details A variety of activities of strongly typed linear algebra with
//! std::mdspan, std::linalg, and mp-units.
[[maybe_unused]] const auto sample{[] {
  // Set up a heterogenous column vector type for the sample.
  using state = column_vector<position, velocity, acceleration>;

  std::vector v0(extents_size<column_extents<3>>, representation{});
  std::mdspan s0{v0.data(), column_extents<3>{}};
  state x0{s0};

  // Elements asignment.
  x0.at<0>(3. * m);
  x0.at<1>(2. * m / s);
  x0.at<2>(1. * m / s2);

  // Printable.
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[3 m], [2 m/s], [1 m/s²]]");

  // Element assignment and access.
  x0.at<1>(2.5 * m / s);
  auto x0_1{x0.at<1>()};
  assert(x0_1 == 2.5 * m / s);
  assert(std::format("{}", x0_1) == "2.5 m/s");

  // Multiplication with a scalar factor.
  scale(3., x0);
  assert(std::format("{}", x0) == "[[9 m], [7.5 m/s], [3 m/s²]]");

  std::vector v4(extents_size<column_extents<3>>, representation{});
  std::mdspan sp4{v4.data(), column_extents<3>{}};
  state x4{sp4};

  // Additions of two vectors of the same types.
  add(x0, x0, x4);
  assert(std::format("{}", x4) == "[[18 m], [15 m/s], [6 m/s²]]");

  using state_transpose = row_vector<position, velocity, acceleration>;

  std::vector v5(extents_size<row_extents<3>>, representation{});
  std::mdspan s5{v5.data(), row_extents<3>{}};
  state_transpose xt5{s5};
  xt5.at<0>(3. * m);
  xt5.at<1>(2. * m / s);
  xt5.at<2>(1. * m / s2);

  using estimate_uncertainty =
      matrix<std::tuple<position, velocity, acceleration>,
             std::tuple<position, velocity, acceleration>>;
  std::vector v6(extents_size<std::extents<std::size_t, 3, 3>>,
                 representation{});
  std::mdspan s6{v6.data(), std::extents<std::size_t, 3, 3>{}};
  estimate_uncertainty p6{s6};

  // Column-row matrix product.
  matrix_product(x0, xt5, p6);
  assert(std::format("{}", p6) == "[[27 m², 18 m²/s, 9 m²/s²],"        //
                                  " [22.5 m²/s, 15 m²/s², 7.5 m²/s³]," //
                                  " [9 m²/s², 6 m²/s³, 3 m²/s⁴]]");

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
