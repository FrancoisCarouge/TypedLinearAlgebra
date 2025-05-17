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

#include "fcarouge/typed_linear_algebra.hpp"

// #include <cassert>
#include <print>

#include <Eigen/Eigen>
#include <mp-units/format.h>
#include <mp-units/framework/quantity.h>
#include <mp-units/framework/quantity_point.h>
#include <mp-units/math.h>
#include <mp-units/systems/isq/thermodynamics.h>
#include <mp-units/systems/si.h>

namespace fcarouge::typed_linear_algebra_internal {

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
};

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
};

//! @brief Multiplies specialization type for uncertainty type deduction.
template <auto Reference>
struct multiplies<mp_units::quantity_point<Reference>, int> {
  [[nodiscard]] inline constexpr auto
  operator()(const mp_units::quantity_point<Reference> &lhs,
             int rhs) const -> mp_units::quantity_point<Reference>;
};

//! @brief Multiplies specialization type for uncertainty type deduction.
template <auto Reference>
struct multiplies<int, mp_units::quantity_point<Reference>> {
  [[nodiscard]] inline constexpr auto
  operator()(int lhs, const mp_units::quantity_point<Reference> &rhs) const
      -> mp_units::quantity_point<Reference>;
};

template <auto Reference1, auto Reference2>
struct multiplies<mp_units::quantity_point<Reference1>,
                  mp_units::quantity_point<Reference2>> {
  [[nodiscard]] inline constexpr auto
  operator()(const mp_units::quantity_point<Reference1> &lhs,
             const mp_units::quantity_point<Reference2> &rhs) const
      -> mp_units::quantity<Reference1 * Reference2>;
};

} // namespace fcarouge::typed_linear_algebra_internal

namespace fcarouge::sample {
namespace {

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

using position = mp_units::quantity_point<mp_units::isq::length[m]>;
using velocity = mp_units::quantity_point<mp_units::isq::velocity[m / s]>;
using acceleration =
    mp_units::quantity_point<mp_units::isq::acceleration[m / s2]>;

template <typename Representation, typename... Types>
using vector =
    typed_column_vector<Eigen::Vector<Representation, sizeof...(Types)>,
                        Types...>;

using state = vector<double, position, velocity, acceleration>;

//! @brief
[[maybe_unused]] auto sample{[] {
  state x0{position{3. * m}, velocity{2. * m / s}, acceleration{1. * m / s2}};

  std::println("X0 = {}", x0);

  state dx{position{3. * m}, velocity{2. * m / s}, acceleration{1. * m / s2}};
  state x1{x0 + dx};

  // a + b - addition where both arguments should be of the same quantity kind
  // and character a - b - subtraction where both arguments should be of the
  // same quantity kind and character a % b - modulo where both arguments should
  // be of the same quantity kind and character a * b - multiplication where one
  // of the arguments has to be a scalar a / b - division where the divisor has
  // to be scalar a ⋅ b - dot product of two vectors a × b - cross product of
  // two vectors |a| - magnitude of a vector a ⊗ b - tensor product of two
  // vectors or tensors a ⋅ b - inner product of two tensors a ⋅ b - inner
  // product of tensor and vector a : b - scalar product of two tensors

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
