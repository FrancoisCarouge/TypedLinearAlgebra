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

//! @brief Strongly typed linear algebra samples.
//!
//! @details A variety of activities of strongly typed linear algebra with Eigen
//! and mp-units.
[[maybe_unused]] auto sample{[] {
  using mp_units::si::unit_symbols::m;
  using mp_units::si::unit_symbols::m2;
  using mp_units::si::unit_symbols::s;
  using mp_units::si::unit_symbols::s2;
  using mp_units::si::unit_symbols::s3;
  constexpr auto s4{pow<4>(s)};
  using mp_units::one;
  using mp_units::si::unit_symbols::A;
  using mp_units::si::unit_symbols::mol;
  using position = quantity<mp_units::isq::length[m]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;
  using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;
  using state = column_vector<position, velocity, acceleration>;

  // Column-vector declaration.
  state x0{3. * m, 2. * m / s, 1. * m / s2};

  // Printable.
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[3 m], [2 m/s], [1 m/s²]]");

  // Element assignment and access.
  x0.at<1>() = 2.5 * m / s;
  auto x0_1{x0.at<1>()};
  assert(x0_1 == 2.5 * m / s);
  assert(std::format("{}", x0_1) == "2.5 m/s");

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

  // Row-vector declaration.
  state_transpose xt5{3. * m, 2. * m / s, 1. * m / s2};
  assert(std::format("{}", xt5) == "[3 m, 2 m/s, 1 m/s²]");

  // Compatible matrix multiplication.
  assert(std::format("{}", x5 * xt5) == "[[9 m², 6 m²/s, 3 m²/s²],"    //
                                        " [6 m²/s, 4 m²/s², 2 m²/s³]," //
                                        " [3 m²/s², 2 m²/s³, 1 m²/s⁴]]");

  // Singleton matrix declaration, for example, but perhaps not a recommended
  // replacement for what should normally just be a `quantity{1. * A / mol}`.
  matrix<std::tuple<decltype(1. * A)>, std::tuple<decltype(1. / mol)>> s1{
      1. * A / mol};
  assert(std::format("{}", s1) == "1 A/mol");

  // Ways to access the singleton matrix.
  s1.at<0, 0>() = 23. * A / mol;
  assert((s1.at<0, 0>() == 23. * A / mol));

  s1(0, 0) = 22. * A / mol;
  assert((s1(0, 0) == 22. * A / mol));

  s1[0, 0] = 21. * A / mol;
  assert((s1[0, 0] == 21. * A / mol));

  s1.at<0>() = 13. * A / mol;
  assert(s1.at<0>() == 13. * A / mol);

  s1(0) = 12. * A / mol;
  assert((s1(0) == 12. * A / mol));

  s1[0] = 11. * A / mol;
  assert((s1[0] == 11. * A / mol));

  s1.at() = 3. * A / mol;
  assert(s1.at() == 3. * A / mol);

  s1() = 2. * A / mol;
  assert(s1() == 2. * A / mol);

  s1[] = 1. * A / mol;
  assert(s1[] == 1. * A / mol);

  // The singleton element can be accessed by conversion to its element type.
  using e = decltype(s1)::element<0, 0>;
  assert(e{s1} == 1. * A / mol);
  // But the element type may not be the same as could be expected.
  static_assert(
      not std::is_same_v<decltype(A / mol), decltype(s1)::element<0, 0>>);

  // More forms of multiplication with a scalar factor.
  assert(std::format("{}", x5 * 2.) == "[[6 m], [4 m/s], [2 m/s²]]");
  assert(std::format("{}", 2. * x5) == "[[6 m], [4 m/s], [2 m/s²]]");
  assert(std::format("{}", xt5 * 2.) == "[6 m, 4 m/s, 2 m/s²]");
  assert(std::format("{}", 2. * xt5) == "[6 m, 4 m/s, 2 m/s²]");
  assert(std::format("{}", (x5 * xt5) * 2.) == "[[18 m², 12 m²/s, 6 m²/s²],"
                                               " [12 m²/s, 8 m²/s², 4 m²/s³],"
                                               " [6 m²/s², 4 m²/s³, 2 m²/s⁴]]");
  assert(std::format("{}", 2. * (x5 * xt5)) == "[[18 m², 12 m²/s, 6 m²/s²],"
                                               " [12 m²/s, 8 m²/s², 4 m²/s³],"
                                               " [6 m²/s², 4 m²/s³, 2 m²/s⁴]]");

  // More forms of multiplication with a strongly typed factor.
  assert(std::format("{}", x5 * (2. * m)) == "[[6 m²], [4 m²/s], [2 m²/s²]]");
  assert(std::format("{}", (2. * m) * x5) == "[[6 m²], [4 m²/s], [2 m²/s²]]");
  assert(std::format("{}", xt5 * (2. * m)) == "[6 m², 4 m²/s, 2 m²/s²]");
  assert(std::format("{}", (2. * m) * xt5) == "[6 m², 4 m²/s, 2 m²/s²]");

  // More forms of multiplication with typed matrices.
  assert(std::format("{}", x5 * s1) == "[[3 A m/mol],"
                                       " [2 A m mol⁻¹ s⁻¹],"
                                       " [1 A m mol⁻¹ s⁻²]]");
  assert(std::format("{}", s1 * xt5) ==
         "[3 A m/mol, 2 A m mol⁻¹ s⁻¹, 1 A m mol⁻¹ s⁻²]");

  // More forms of division with a scalar term.
  assert(std::format("{}", x5 / 2.) == "[[1.5 m], [1 m/s], [0.5 m/s²]]");
  assert(std::format("{}", 1. / x5) ==
         "[0.3333333333333333 1/m, 0 s/m, 0 s²/m]");
  assert(std::format("{}", xt5 / 2.) == "[1.5 m, 1 m/s, 0.5 m/s²]");

  // More forms of division with a strongly typed factor.
  assert(std::format("{}", x5 / (2. * m)) == "[[1.5], [1 1/s], [0.5 1/s²]]");
  assert(std::format("{}", (1. * m) / x5) == "[0.3333333333333333, 0 s, 0 s²]");
  assert(std::format("{}", xt5 / (2. * m)) == "[1.5, 1 1/s, 0.5 1/s²]");

  // More forms of division with typed matrices.
  assert(std::format("{}", s1 / s1) == "1");
  assert(std::format("{}", x5 / x5) == "[[1, 0 s, 0 s²],"
                                       " [0.6666666666666666 1/s, 0, 0 s],"
                                       " [0.3333333333333333 1/s², 0 1/s, 0]]");

  //! @todo A few more expression may be valid, mising: x5 / s1? s1 / x5?

  // Homogeneously quantity typed vector.
  using vector3d =
      column_vector<representation, representation, representation>;
  using velocity3d = column_vector<velocity, velocity, velocity>;

  vector3d v{1., 2., 3.};
  velocity3d v0{v * m / s};
  assert(std::format("{}", v0) == "[[1 m/s], [2 m/s], [3 m/s]]");

  velocity a[]{1. * m / s, 2. * m / s, 3. * m / s};
  velocity3d v1{a};
  assert(std::format("{}", v1) == "[[1 m/s], [2 m/s], [3 m/s]]");

  // Vector typed access.
  v0.at<1>() = 4. * m / s;
  assert(v0.at<1>() == 4. * m / s);

  // Vector and uniform typed access.
  v0[1] = 3. * m / s;
  assert(v0[1] == 3. * m / s);

  v0(1) = 2. * m / s;
  assert(v0(1) == 2. * m / s);

  // Beware of non-evaluated template expression: these types are not the same.
  auto a0{vector3d{1., 2., 3.} * mp_units::isq::velocity[m / s]};
  static_assert(not std::is_same_v<decltype(a0), velocity3d>);

  //! @todo Mp-units seems to have support for vector quantity? Of the form:
  //! mp_units::Quantity auto v0{...};

  // Addition where both arguments should be of the same quantity kind and
  // character.
  //! @todo Add additional valid forms for other matrices.
  assert(std::format("{}", v0 + v0) == "[[2 m/s], [4 m/s], [6 m/s]]");

  // Subtraction where both arguments should be of the same quantity kind and
  // character.
  //! @todo Add additional valid forms for other matrices.
  assert(std::format("{}", v0 - v0) == "[[0 m/s], [0 m/s], [0 m/s]]");

  // Matrix and uniform typed access.
  using position_2d_uncertainty =
      matrix<std::tuple<position, position>, std::tuple<position, position>>;
  position_2d_uncertainty p0;

  p0[0, 1] = 9. * m2;
  assert((p0[0, 1] == 9. * m2));

  p0(0, 1) = 16. * m2;
  assert((p0(0, 1) == 16. * m2));

  //! @todo Modulo where both arguments should be of the same quantity kind and
  //! character.
  //! @todo Dot product of two vectors: a ⋅ b.
  //! @todo Cross product of two vectors: a × b.
  //! @todo Magnitude of a vector: |a|.
  //! @todo Tensor product of two vectors or tensors: a ⊗ b.
  //! @todo Inner product of two tensors: a ⋅ b.
  //! @todo Inner product of tensor and vector: a ⋅ b.
  //! @todo Scalar product of two tensors: a : b.
  //! @todo Further compatibility with STL/Ranges/std::linalg: sum_v=accumulate(
  //! v.cbegin()+1,v.cend()-1,0.0)
  //! @todo prod_v=accumulate(v.cbegin(),v.cend(),1.0,std::multiplies<float>())
  //! @todo is_all_zero=std::ranges::all_of(v,[](floatf){returnf==0;})
  //! @todo count_v=std::ranges::count_if(v,[](floatf){returnf>5.0;})
  //! @todo max_v=*std::ranges::max_element(v)
  //! @todo std::ranges::fill(v,5.0)
  //! @todo std::ranges::copy(v,w.begin())
  //! @todo std::ranges::copy_if(v,w.begin(),[](floatf){returnf>5.0;})
  //! @todo std::ranges::transform(v,w.begin(),[](floatf){returnf+1.0;})
  //! @todo x=inner_product(v.cbegin(),v.cend(),w.begin(),0.0)

  // 1-D vehicle location Kalman estimation.
  state x{0. * m, 0. * m / s, 0. * m / s2};
  std::println("X: {}", x);
  // X: [[0 m],
  //     [0 m/s],
  //     [0 m/s²]]

  using estimate_uncertainty =
      matrix<std::tuple<position, velocity, acceleration>,
             std::tuple<position, velocity, acceleration>>;
  estimate_uncertainty p;
  p.at<0, 0>() = 500. * m2;
  p.at<1, 1>() = 500. * m2 / s2;
  p.at<2, 2>() = 500. * m2 / s4;
  std::println("P: {}", p);
  // P: [[500 m²,     0 m²/s,    0 m²/s²],
  //     [  0 m²/s, 500 m²/s²,   0 m²/s³],
  //     [  0 m²/s²,  0 m²/s³, 500 m²/s⁴]]

  using process_uncertainty = estimate_uncertainty;
  process_uncertainty q;
  q.at<0, 0>() = 0.01 * m2;
  q.at<0, 1>() = 0.02 * m2 / s;
  q.at<0, 2>() = 0.02 * m2 / s2;
  q.at<1, 0>() = 0.02 * m2 / s;
  q.at<1, 1>() = 0.04 * m2 / s2;
  q.at<1, 2>() = 0.04 * m2 / s3;
  q.at<2, 0>() = 0.02 * m2 / s2;
  q.at<2, 1>() = 0.04 * m2 / s3;
  q.at<2, 2>() = 0.04 * m2 / s4;
  std::println("Q: {}", q);
  // Q: [[0.01 m²,    0.02 m²/s,  0.02 m²/s²],
  //     [0.02 m²/s,  0.04 m²/s², 0.04 m²/s³],
  //     [0.02 m²/s², 0.04 m²/s³, 0.04 m²/s⁴]]

  using output_uncertainty = quantity<m2>;
  output_uncertainty r{9. * m2};
  std::println("R: {}", r);
  // R: 9 m²

  using output_model = row_vector<quantity<one>, quantity<s>, quantity<s2>>;
  output_model h;
  h.at<0, 0>() = 1.;
  std::println("H: {}", h);
  // H: [1, 0 s, 0 s²]

  using state_transition =
      matrix<std::tuple<position, velocity, acceleration>,
             std::tuple<quantity<one / m>, quantity<s / m>, quantity<s2 / m>>>;
  state_transition f;
  f.at<0, 0>() = 1.;
  f.at<0, 1>() = 1. * s;
  f.at<0, 2>() = 0.5 * s2;
  f.at<1, 1>() = 1.;
  f.at<1, 2>() = 1. * s;
  f.at<2, 2>() = 1.;
  std::println("F: {}", f);
  // F: [[1, 1 s, 0.5 s²],
  //     [0 1/s, 1, 1 s],
  //     [0 1/s², 0 1/s, 1]]

  // Prediction stage of the filter estimated state.
  x = f * x;
  p = f * p * transposed(f) + q;

  // Update stage of the filter from output measurements.
  using output = position;
  output z{-393.66 * m};

  using innovation_uncertainty = output_uncertainty;
  innovation_uncertainty si{h * p * transposed(h) + r};

  using unevaluated_gain =
      decltype(std::declval<state>() / std::declval<output>());
  using gain =
      matrix<unevaluated_gain::row_indexes, unevaluated_gain::column_indexes>;
  gain k{p * transposed(h) / si};

  using innovation = output;
  innovation y{z - h * x};
  x = x + k * y;

  std::println("X: {}", x);
  // X: [[-390.53 m],
  //     [-260.36 m/s],
  //     [ -86.79 m/s²]]

  using unevaluated_kh =
      decltype(std::declval<gain>() * std::declval<output_model>());
  using kh =
      matrix<unevaluated_kh::row_indexes, unevaluated_kh::column_indexes>;
  kh i;
  i.at<0, 0>() = 1.;
  i.at<1, 1>() = 1.;
  i.at<2, 2>() = 1.;
  p = (i - k * h) * p * transposed(i - k * h) + k * r * transposed(k);
  std::println("P: {}", p);
  // P: [[8.92 m²,      5.95 m²/s,    1.98 m²/s²],
  //     [5.95 m²/s,  503.98 m²/s², 334.67 m²/s³],
  //     [1.98 m²/s², 334.67 m²/s³, 444.91 m²/s⁴]]

  return 0;
}()};
} // namespace
} // namespace sample
} // namespace fcarouge
