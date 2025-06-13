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

namespace fcarouge {
template <mp_units::Quantity To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] inline constexpr To operator()(const From &value) const {
    return value * To::reference;
  }
};

template <typename To, mp_units::Quantity From>
struct element_caster<To, From> {
  [[nodiscard]] inline constexpr To operator()(const From &value) const {
    return value.numerical_value_in(value.unit);
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To &, From &> {
  [[nodiscard]] inline constexpr To &operator()(From &value) const {
    return reinterpret_cast<To &>(value);
  }
};

// Make the useful tools available to the user. Revise evaluate, internalize in
// the matrix?

// template <mp_units::Quantity Denominator, typename Matrix, typename
// RowIndexes,
//           typename ColumnIndexes>
// [[nodiscard]] inline constexpr auto
// operator/(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
//           const Denominator& rhs) {
//   return typed_matrix<tla::evaluate<Matrix>, RowIndexes, ColumnIndexes>{
//       lhs.data() / rhs};
// }

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
  using mp_units::si::unit_symbols::m;
  using mp_units::si::unit_symbols::s;
  using mp_units::si::unit_symbols::s2;

  using position = quantity<mp_units::isq::length[m]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;
  using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

  using state = column_vector<position, velocity, acceleration>;

  // Declaration.
  state x0{3. * m, 2. * m / s, 1. * m / s2};
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[3 m], [2 m/s], [1 m/s²]]");

  // Element access.
  x0.at<1>() = 2.5 * m / s;
  assert(x0.at<1>() == 2.5 * m / s);
  assert(std::format("{}", x0.at<1>()) == "2.5 m/s");

  // Multiplication with a scalar.
  state x1{x0 * 3.};
  std::println("x1 = {}", x1);
  assert(std::format("{}", x1) == "[[9 m], [7.5 m/s], [3 m/s²]]");

  // Division with a scalar divisor.
  state x2{x1 / 2.};
  std::println("x2 = {}", x2);
  assert(std::format("{}", x2) == "[[4.5 m], [3.75 m/s], [1.5 m/s²]]");

  // Substraction of two vectors of the same types.
  state x3{x2 - x0};
  std::println("x3 = {}", x3);
  assert(std::format("{}", x3) == "[[1.5 m], [1.25 m/s], [0.5 m/s²]]");

  // Substraction of two vectors of the same types.
  state x4{x3 + x3};
  std::println("x4 = {}", x4);
  assert(std::format("{}", x4) == "[[3 m], [2.5 m/s], [1 m/s²]]");

  using state_transpose = row_vector<position, velocity, acceleration>;

  state x5{3. * m, 2. * m / s, 1. * m / s2};
  state_transpose xt5{3. * m, 2. * m / s, 1. * m / s2};
  using uncertainty = matrix<std::tuple<position, velocity, acceleration>,
                             std::tuple<position, velocity, acceleration>>;

  // Inner product.
  double d5{xt5 * x5};
  std::println("d5 = {}", d5);
  assert(std::format("{}", d5) == "14");

  using uncertainty = matrix<std::tuple<position, velocity, acceleration>,
                             std::tuple<position, velocity, acceleration>>;

  // Column-vector multiplication.
  uncertainty p5{x5 * xt5};
  std::println("p5 = {}", p5);
  assert(std::format("{}", p5) == "[[9 m², 6 m²/s, 3 m²/s²],"    //
                                  " [6 m²/s, 4 m²/s², 2 m²/s³]," //
                                  " [3 m²/s², 2 m²/s³, 1 m²/s⁴]]");

  std::println("p6 = {}", p5 * p5); // Is this correct? It can't be.

  std::println("p7 = {}", p5 * 2);
  std::println("p8 = {}", 2 * p5);
  // std::println("p7 = {}", p5 * (2. * m));

  // TODO Continue sample with temperatures to show the gaps with quantity
  // points and show how complicated things get.

  // TODO Compatibility with STL/Ranges.
  // sum_v = accumulate (v. cbegin ()+1 , v. cend () -1 , 0.0);
  // prod_v = accumulate (v . cbegin () , v. cend () , 1.0 , std ::
  // multiplies < float >()); is_all_zero = std :: ranges :: all_of (v , [](
  // float f ){ return f ==0;}); count_v = std :: ranges :: count_if (v ,
  // []( float f ){ return f >5.0;}); max_v = * std :: ranges :: max_element
  // (v ); std :: ranges :: fill (v , 5.0); std :: ranges :: copy (v , w.
  // begin ()); std :: ranges :: copy_if (v , w. begin () , []( float f ){
  // return f >5.0;}); std
  // :: ranges :: transform (v , w. begin () , []( float f ){ return f
  // +1.0;}); x = inner_product (v. cbegin () , v. cend () , w. begin () ,
  // 0.0);

  // a + b - addition where both arguments should be of the same quantity
  // kind and character a - b - subtraction where both arguments should be
  // of the same quantity kind and character a % b - modulo where both
  // arguments should be of the same quantity kind and character a * b -
  // multiplication where one of the arguments has to be a scalar a / b -
  // division where the divisor has to be scalar a ⋅ b - dot product of two
  // vectors a × b - cross product of two vectors |a| - magnitude of a
  // vector a ⊗ b - tensor product of two vectors or tensors a ⋅ b - inner
  // product of two tensors a ⋅ b - inner product of tensor and vector a : b
  // - scalar product of two tensors

  return 0;
}()};
} // namespace
} // namespace sample
} // namespace fcarouge
