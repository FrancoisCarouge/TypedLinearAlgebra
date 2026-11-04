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

#include "fcarouge/linalg.hpp"

#include <cassert>
#include <concepts>
#include <tuple>
#include <type_traits>

namespace fcarouge::test {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::km;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;
using mp_units::si::unit_symbols::m4;
using mp_units::si::unit_symbols::s;

namespace {
//! @test Verifies the sparse storage with mp-units quantities, from the
//! bilinear quantity matrix to the linear map of a state transition.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using length = quantity<mp_units::isq::length[m]>;
  using area = quantity<mp_units::isq::area[m2]>;
  using speed = quantity<mp_units::isq::speed[m / s]>;
  using indexes = std::tuple<length, length>;
  using result_indexes = std::tuple<area, area>;

  using diagonal_matrix =
      sparse_matrix<representation, sparse::diagonal<2>, indexes, indexes>;

  // A read of a structural zero has the element type.
  const diagonal_matrix lhs{sparse::make_matrix<diagonal_matrix>(
      1. * m2, sparse::structural_zero, sparse::structural_zero, 2. * m2)};

  static_assert(std::same_as<std::remove_cvref_t<decltype(lhs.at<0, 1>())>,
                             std::remove_cvref_t<decltype(lhs.at<0, 0>())>>);
  assert((lhs.at<0, 0>() == 1. * m2));
  assert((lhs.at<0, 1>() == 0. * m2));

  // A product computes the element types, a sparse pattern stays sparse.
  const sparse_matrix<representation, sparse::tridiagonal<2>, indexes, indexes>
      rhs{sparse::make_matrix<sparse_matrix<
          representation, sparse::tridiagonal<2>, indexes, indexes>>(
          3. * m2, 4. * m2, 5. * m2, 6. * m2)};
  const sparse_matrix<representation, sparse::dense<2, 2>, result_indexes,
                      result_indexes>
      product{lhs * rhs};

  assert((product.at<0, 0>() == 3. * m4));
  assert((product.at<1, 1>() == 12. * m4));

  const auto squared_value{lhs * lhs};

  static_assert(decltype(squared_value)::matrix::structure.count == 2);
  assert((squared_value.at<1, 1>() == 4. * m4));
  assert((squared_value.at<1, 0>() == 0. * m4));

  // A sum, a difference, a negation, a scaling, a transposition, a norm.
  assert(((lhs + rhs).at<0, 1>() == 4. * m2));
  assert(((lhs - rhs).at<1, 0>() == -5. * m2));
  assert(((-lhs).at<1, 1>() == -2. * m2));
  assert(((lhs * 2.).at<1, 1>() == 4. * m2));
  assert((transposed(rhs).at<0, 1>() == 5. * m2));
  assert(matrix_one_norm(lhs) == 2. * m2);
  assert(matrix_inf_norm(rhs) == 11. * m2);

  // A given unit converts to the element's on construction and on write.
  using position_vector =
      sparse_column_vector<representation, sparse::make_pattern<2, 1>({{0, 0}}),
                           length, length>;

  position_vector position{
      sparse::make_matrix<position_vector>(1. * km, sparse::structural_zero)};

  assert(position.at<0>() == 1000. * m);
  assert(position.at<1>() == 0. * m);

  position.at<0>(2. * km);

  assert(position.at<0>() == 2000. * m);

  // Heterogeneous elements with a structural zero of its own quantity.
  using state =
      sparse_column_vector<representation, sparse::make_pattern<2, 1>({{0, 0}}),
                           length, speed>;

  const state rest{sparse::make_matrix<state>(3. * m, sparse::structural_zero)};

  assert(rest.at<0>() == 3. * m);
  assert(rest.at<1>() == 0. * m / s);

  // A vector quantity over a sparse representation scales its stored elements.
  using vector3d =
      sparse_matrix<double, sparse::make_pattern<3, 1>({{0, 0}, {2, 0}})>;

  const mp_units::quantity velocity{mp_units::isq::velocity(
      sparse::make_matrix<vector3d>(30., sparse::structural_zero, 40.) * m /
      s)};

  assert(magnitude(velocity.numerical_value_in(m / s)) == 50.);
  assert(velocity.in(km / s).numerical_value_in(km / s).at<2>() == 0.04);

  // A state transition, a linear map: the columns are reciprocal to the
  // state's. An upper triangular transition of a position and a speed.
  using inverse_length = quantity<mp_units::one / mp_units::isq::length[m]>;
  using inverse_speed = quantity<mp_units::one / mp_units::isq::speed[m / s]>;

  using transition = sparse_matrix<representation, sparse::upper_triangular<2>,
                                   std::tuple<length, speed>,
                                   std::tuple<inverse_length, inverse_speed>>;

  transition model;

  model.at<0, 0>(1. * mp_units::one);
  model.at<0, 1>(2. * s);
  model.at<1, 1>(1. * mp_units::one);

  assert((model.at<1, 0>() == 0. / s));

  const sparse_column_vector<representation, sparse::dense<2, 1>, length, speed>
      initial{1. * m, 3. * m / s};
  const sparse_column_vector<representation, sparse::dense<2, 1>, length, speed>
      final_state{model * initial};

  assert(final_state.at<0>() == 7. * m);
  assert(final_state.at<1>() == 3. * m / s);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
