/* Typed Linear Algebra
Version 0.1.0
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

#include "fcarouge/typed_linear_algebra.hpp"

#include <cstddef>
#include <format>
#include <linalg>
#include <mdspan>
#include <print>
#include <tuple>
#include <type_traits>
#include <vector>

#include <mp-units/framework/quantity.h>
#include <mp-units/framework/quantity_point.h>
#include <mp-units/math.h>
#include <mp-units/systems/isq/thermodynamics.h>
#include <mp-units/systems/si.h>

namespace fcarouge {
// Teach the typed linear algebra library how to convert Eigen' underlying
// scalar types to and from mp-units' types.
template <typename To, mp_units::Quantity From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(From value) const {
    return value.numerical_value_in(value.unit);
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(From value) const {
    return value * To::reference;
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To &, From &> {
  [[nodiscard]] constexpr To &operator()(From &value) const {
    return reinterpret_cast<To &>(value);
  }
};

template <typename To, mp_units::Reference From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(From) const { return 1.; }
};

namespace sample {
namespace {
// Set up heterogenously unit typed linear algebra types.
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

template <mp_units::Reference auto QuantityReference>
using quantity_point =
    mp_units::quantity_point<QuantityReference,
                             mp_units::default_point_origin(QuantityReference),
                             representation>;

template <typename RowIndexes, typename ColumnIndexes>
using matrix = typed_matrix<
    std::mdspan<representation,
                std::extents<std::size_t, std::tuple_size_v<RowIndexes>,
                             std::tuple_size_v<ColumnIndexes>>>,
    RowIndexes, ColumnIndexes>;

template <typename... Types>
using column_vector = typed_column_vector<
    std::mdspan<representation, std::extents<std::size_t, sizeof...(Types), 1>>,
    Types...>;

template <typename... Types>
using row_vector = typed_row_vector<
    std::mdspan<representation, std::extents<std::size_t, 1, sizeof...(Types)>>,
    Types...>;

template <std::size_t Rows>
using column_extents = std::extents<std::size_t, Rows, 1>;

template <std::size_t Columns>
using row_extents = std::extents<std::size_t, 1, Columns>;

template <typename Extents>
constexpr std::size_t extents_size{[] {
  std::size_t size{1};
  Extents extents{};
  for (std::size_t i{0}; i < extents.rank(); ++i) {
    size *= extents.extent(i);
  }
  return size;
}()};

//! @brief Strongly typed linear algebra samples.
//!
//! @details A variety of activities of strongly typed linear algebra with
//! std::mdspan, std::linalg, and mp-units.
[[maybe_unused]] auto sample{[] {
  using mp_units::si::unit_symbols::m;
  using mp_units::si::unit_symbols::m2;
  using mp_units::si::unit_symbols::s;
  using mp_units::si::unit_symbols::s2;
  using mp_units::si::unit_symbols::s3;
  using mp_units::one;
  using mp_units::si::unit_symbols::A;
  using mp_units::si::unit_symbols::mol;
  using position = quantity<mp_units::isq::length[m]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;
  using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;
  using state = column_vector<position, velocity, acceleration>;

  std::vector v0(extents_size<column_extents<3>>, representation{});
  std::mdspan s0{v0.data(), column_extents<3>{}};
  state x0{s0};

  // Elements asignment.
  x0.at<0>() = 3. * m;
  x0.at<1>() = 2. * m / s;
  x0.at<2>() = 1. * m / s2;

  // Printable.
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[3 m], [2 m/s], [1 m/s²]]");

  // Element assignment and access.
  x0.at<1>() = 2.5 * m / s;
  auto x0_1{x0.at<1>()};
  assert(x0_1 == 2.5 * m / s);
  assert(std::format("{}", x0_1) == "2.5 m/s");

  // Multiplication with a scalar factor.
  scale(3., x0);
  assert(std::format("{}", x0) == "[[9 m], [7.5 m/s], [3 m/s²]]");

  std::vector v4(extents_size<column_extents<3>>, representation{});
  std::mdspan s4{v4.data(), column_extents<3>{}};
  state x4{s4};

  // Additions of two vectors of the same types.
  add(x0, x0, x4);
  assert(std::format("{}", x4) == "[[18 m], [15 m/s], [6 m/s²]]");

  using state_transpose = row_vector<position, velocity, acceleration>;

  std::vector v5(extents_size<row_extents<3>>, representation{});
  std::mdspan s5{v5.data(), row_extents<3>{}};
  state_transpose xt5{s5};
  xt5.at<0>() = 3. * m;
  xt5.at<1>() = 2. * m / s;
  xt5.at<2>() = 1. * m / s2;

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
} // namespace sample
} // namespace fcarouge
