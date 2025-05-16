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

#include "fcarouge/linalg.hpp"
#include "fcarouge/typed_linear_algebra.hpp"

#include <format>
#include <print>
#include <tuple>
#include <type_traits>

#include <Eigen/Eigen>
#include <mp-units/framework/quantity.h>
#include <mp-units/framework/quantity_point.h>
#include <mp-units/math.h>
#include <mp-units/systems/isq/thermodynamics.h>
#include <mp-units/systems/si.h>

namespace fcarouge {
// Teach the typed linear algebra library about mp-units type conversions.
template <mp_units::Quantity To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(const From &value) const {
    return value * To::reference;
  }
};

template <typename To, mp_units::Quantity From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(const From &value) const {
    return value.numerical_value_in(value.unit);
  }
};

template <typename To, mp_units::Reference From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(const From &) const { return 1.; }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To &, From &> {
  [[nodiscard]] constexpr To &operator()(From &value) const {
    return reinterpret_cast<To &>(value);
  }
};

namespace sample {
namespace {
// Set up heterogenously typed linear algebra types.
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

  // Printable.
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[3 m], [2 m/s], [1 m/s²]]");

  // Element access.
  x0.at<1>() = 2.5 * m / s;
  assert(x0.at<1>() == 2.5 * m / s);
  assert(std::format("{}", x0.at<1>()) == "2.5 m/s");

  // Multiplication with a scalar factor.
  state x1{x0 * 3.};
  assert(std::format("{}", x1) == "[[9 m], [7.5 m/s], [3 m/s²]]");

  // Division with a scalar divisor.
  state x2{x1 / 2.};
  assert(std::format("{}", x2) == "[[4.5 m], [3.75 m/s], [1.5 m/s²]]");

  // Substraction of two vectors of the same types.
  state x3{x2 - x0};
  assert(std::format("{}", x3) == "[[1.5 m], [1.25 m/s], [0.5 m/s²]]");

  // Additions of two vectors of the same types.
  state x4{x3 + x3};
  assert(std::format("{}", x4) == "[[3 m], [2.5 m/s], [1 m/s²]]");

  state x5{3. * m, 2. * m / s, 1. * m / s2};

  // Multiplication with a strongly typed factor.
  assert(std::format("{}", x5 * (2. * m)) == "[[6 m²], [4 m²/s], [2 m²/s²]]");
  assert(std::format("{}", (0.5 / m) * x5) == "[[1.5], [1 1/s], [0.5 1/s²]]");

  using state_transpose = row_vector<position, velocity, acceleration>;
  state_transpose xt5{3. * m, 2. * m / s, 1. * m / s2};

  // Matrix multiplication.
  assert(std::format("{}", x5 * xt5) == "[[9 m², 6 m²/s, 3 m²/s²],"    //
                                        " [6 m²/s, 4 m²/s², 2 m²/s³]," //
                                        " [3 m²/s², 2 m²/s³, 1 m²/s⁴]]");

  // Tests:
  std::println("====== TEST =======");
  // column_vector<position> s1{1. * m};
  using mp_units::si::unit_symbols::A;
  using mp_units::si::unit_symbols::mol;
  matrix<std::tuple<decltype(1. * A)>, std::tuple<decltype(1. / mol)>> s1{
      1. * A / mol};
  std::println("00: S1: \n    {}", s1);
  std::println("01: X5: \n    {}", x5);
  // [[3 m],
  //  [2 m/s],
  //  [1 m/s²]]

  std::println("01b: XT5: \n    {}", xt5);

  std::println("02: X5 * 1.: \n    {}", x5 * 1.);
  std::println("03: 1. * X5: \n    {}", 1. * x5);
  std::println("04: XT5 * 1.: \n    {}", xt5 * 1.);
  std::println("05: 1. * XT5: \n    {}", 1. * xt5);

  std::println("06: X5 * 1. m: \n    {}", x5 * (1. * m));
  std::println("07: 1. m * X5: \n    {}", (1. * m) * x5);
  std::println("08: XT5 * 1. m: \n    {}", xt5 * (1. * m));
  std::println("09: 1. m * XT5: \n    {}", (1. * m) * xt5);

  std::println("10: X5 / 1.: \n    {}", x5 / 1.);
  std::println("11: 1. / X5: \n    {}", 1. / x5);
  std::println("12: XT5 / 1.: \n    {}", xt5 / 1.);

  std::println("14: X5 / 1. m: \n    {}", x5 / (1. * m));
  std::println("15: 1. m / X5: \n    {}", (1. * m) / x5);
  std::println("16: XT5 / 1. m: \n    {}", xt5 / (1. * m));

  std::println("17a: s1 / s1: \n    {}", s1 / s1);
  // std::println("17b: s1 / X5: \n    {}", s1 / x5); // Non-sense?
  // std::println("17c: XT5 / s1: \n    {}", xt5 / s1);

  std::println("19: X5 x S1: \n    {}", x5 * s1);
  std::println("20: S1 x XT5: \n    {}", s1 * xt5);
  std::println("21: X5 x XT5: \n    {}", x5 * xt5);

  // [[9 m²,    6 m²/s,  3 m²/s²],
  //  [6 m²/s,  4 m²/s², 2 m²/s³],
  //  [3 m²/s², 2 m²/s³, 1 m²/s⁴]]

  std::println("22: (X5 x S1) x (S1 x XT5): \n    {}", (x5 * s1) * (s1 * xt5));

  std::println("23: X5 / X5: \n    {}", x5 / x5);

  // EXPECTED?
  //                         1        s     s2
  // 1    [[                 1,       0 s,  0 s²],
  // 1/s  [0.6666666666666666 1/s,   0,    0 s],
  // 1/s2 [0.3333333333333333 1/s²,  0 1/s, 0]]

  std::println("24: (X5 / X5) * X5: \n    {}", (x5 / x5) * x5);

  // error: static assertion failed due to requirement '
  // std::is_same_v<
  // mp_units::quantity<mp_units::reference<
  //   mp_units::isq::velocity, mp_units::derived_unit<mp_units::si::metre,
  //   mp_units::per<mp_units::si::second>>>{}, double>,
  // mp_units::quantity<mp_units::reference<
  //   mp_units::isq::acceleration, mp_units::derived_unit<mp_units::si::metre,
  //   mp_units::per<mp_units::power<mp_units::si::second, 2>>>>{}, double>>
  // ': Matrix multiplication requires compatible types.

  //         1
  // m    [[3 m],
  // m/s  [2 m/s],
  // m/s2 [1 m/s²]]

  // std::println("23: (X5 x XT5) / XT5 ~= X5: \n    {}",
  //              (x5 * xt5) / xt5); // WRONG

  // /
  //  [3 m,     2 m/s,   1 m/s²]

  //  76 |     p = estimate_uncertainty{(i - k * h) * p * t(i - k * h) + k * r *
  //  t(k)};
  //     |                                          ^

  std::println("25: \n (k * h) * p   {}", ((x5 / x5) * (x5 / x5)) * (x5 * xt5));

  using vector3d =
      column_vector<representation, representation, representation>;

  // TODO Can we make this also a true mp-units quantity?
  /*mp_units::Quantity */ auto v0{vector3d{1., 2., 3.} *
                                  mp_units::isq::velocity[m / s]};
  assert(std::format("{}", v0) == "[[1 m/s], [2 m/s], [3 m/s]]");

  // Beware of expression templates: these types are not the same.
  using velocity3d = column_vector<velocity, velocity, velocity>;
  static_assert(not std::is_same_v<decltype(v0), velocity3d>);

  // TODO linalg need.
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

  return 0;
}()};
} // namespace
} // namespace sample
} // namespace fcarouge
