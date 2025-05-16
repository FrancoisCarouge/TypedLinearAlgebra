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
//! and Eigen. This library composes Eigen as the linear algebra backend with
//! index typed as mp-units types. This sample explicitly uses double precision
//! floating point numbers. This sample uses Eigen linear algebra as the linear
//! algebra backend. This sample uses mp-units types for the typed linear
//! algebra.

#include "fcarouge/typed_linear_algebra.hpp"

#include <format>
#include <print>
#include <tuple>
#include <type_traits>

#include <Eigen/Eigen>
#include <mp-units/format.h>
#include <mp-units/framework/quantity.h>
#include <mp-units/framework/quantity_point.h>
#include <mp-units/math.h>
#include <mp-units/systems/isq/thermodynamics.h>
#include <mp-units/systems/si.h>

// template <mp_units::QuantityPoint QuantityPoint, typename Char>
// struct std::formatter<QuantityPoint, Char>
//     : formatter<typename QuantityPoint::quantity_type, Char> {
//   template <typename FormatContext>
//   constexpr auto format(const QuantityPoint &value,
//                         FormatContext &format_context) const
//       -> FormatContext::iterator {
//     formatter<typename QuantityPoint::quantity_type, Char>::format(
//         value - value.absolute_point_origin, format_context);
//     return format_context.out();
//   }
// };

namespace fcarouge {
namespace typed_linear_algebra_internal {
//! @todo Users must not have to specialize internal types for this library to
//! work. Find a solution.
template <auto Reference, auto OriginPoint, typename Representation>
struct element_traits<
    Representation,
    mp_units::quantity_point<Reference, OriginPoint, Representation>> {
  [[nodiscard]] static inline constexpr Representation to_underlying(
      const mp_units::quantity_point<Reference, OriginPoint, Representation>
          &value) {
    return value.quantity_from(OriginPoint).numerical_value_in(value.unit);
  }

  [[nodiscard]] static inline constexpr mp_units::quantity_point<
      Reference, OriginPoint, Representation> &
  from_underlying(Representation &value) {
    return reinterpret_cast<
        mp_units::quantity_point<Reference, OriginPoint, Representation> &>(
        value);
  }

  [[nodiscard]] static inline constexpr mp_units::quantity_point<
      Reference, OriginPoint, Representation>
  from_underlying(const Representation &value) {
    return mp_units::quantity_point<Reference, OriginPoint, Representation>{
        value * Reference};
  }
};

// Can we do without point in this sample?
template <auto Reference, typename Representation>
struct element_traits<Representation,
                      mp_units::quantity<Reference, Representation>> {
  [[nodiscard]] static inline constexpr Representation
  to_underlying(const mp_units::quantity<Reference, Representation> &value) {
    return value.numerical_value_in(value.unit);
  }

  [[nodiscard]] static inline constexpr mp_units::quantity<Reference,
                                                           Representation> &
  from_underlying(Representation &value) {
    return reinterpret_cast<mp_units::quantity<Reference, Representation> &>(
        value);
  }

  [[nodiscard]] static inline constexpr mp_units::quantity<Reference,
                                                           Representation>
  from_underlying(const Representation &value) {
    return mp_units::quantity<Reference, Representation>{value * Reference};
  }
};

} // namespace typed_linear_algebra_internal
namespace sample {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

template <mp_units::Reference auto QuantityReference>
using quantity_point =
    mp_units::quantity_point<QuantityReference,
                             mp_units::default_point_origin(QuantityReference),
                             representation>;

template <typename RowIndexes, typename ColumnIndexes>
using matrix =
    typed_matrix<Eigen::Matrix<representation, std::tuple_size_v<RowIndexes>,
                               std::tuple_size_v<ColumnIndexes>>,
                 RowIndexes, ColumnIndexes>;

template <typename... Types>
using column_vector =
    typed_column_vector<Eigen::Vector<representation, sizeof...(Types)>,
                        Types...>;

template <typename... Types>
using row_vector =
    typed_row_vector<Eigen::RowVector<representation, sizeof...(Types)>,
                     Types...>;

[[maybe_unused]] auto sample{[] {
  using mp_units::delta;
  using mp_units::si::unit_symbols::m;
  using mp_units::si::unit_symbols::s;
  using mp_units::si::unit_symbols::s2;

  using position = quantity<mp_units::isq::length[mp_units::si::metre]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;
  using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

  using state = column_vector<position, velocity, acceleration>;

  state x0{3. * m, 2. * m / s, 1. * m / s2};
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[3 m], [2 m/s], [1 m/s²]]");

  x0.at<1>() = 2.5 * m / s;
  assert(x0.at<1>() == 2.5 * m / s);
  assert(std::format("{}", x0.at<1>()) == "2.5 m/s");

  state x1{x0 * 3.};
  std::println("x1 = {}", x1);
  assert(std::format("{}", x1) == "[[9 m], [7.5 m/s], [3 m/s²]]");

  state x2{x1 / 2.};
  std::println("x2 = {}", x2);
  assert(std::format("{}", x2) == "[[4.5 m], [3.75 m/s], [1.5 m/s²]]");

  state x3{x2 - x0};
  std::println("x3 = {}", x3);
  assert(std::format("{}", x3) == "[[1.5 m], [1.25 m/s], [0.5 m/s²]]");

  state x4{x3 + x3};
  std::println("x4 = {}", x4);
  assert(std::format("{}", x4) == "[[3 m], [2.5 m/s], [1 m/s²]]");

  // state x5{x4 * (1. * m)};
  // std::println("x5 = {}", x5);
  // assert(std::format("{}", x5) == "[[3 m], [2.5 m/s], [1 m/s²]]");

  // TODO Continue sample with temperatures to show the gaps with quantity
  // points and show how complicated things get.

  return 0;
}()};
} // namespace
} // namespace sample
} // namespace fcarouge
