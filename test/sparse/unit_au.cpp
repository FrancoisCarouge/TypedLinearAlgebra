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

#include <au/units/meters.hh>
#include <au/units/seconds.hh>

#include <cassert>
#include <concepts>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the sparse storage with Au quantities: the per-element typing
//! of the indexes, whatever the pattern, a structural zero included.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using au::symbols::m;
  using au::symbols::s;

  constexpr auto m2{au::squared(m)};
  constexpr auto m4{pow<4>(m)};

  using length = au::QuantityD<au::Meters>;
  using duration_type = au::QuantityD<au::Seconds>;
  using area = au::QuantityD<au::UnitPowerT<au::Meters, 2>>;
  using indexes = std::tuple<length, length>;
  using result_indexes = std::tuple<area, area>;

  // A uniform quantity matrix, a read of a structural zero has the element
  // type.
  using diagonal_matrix =
      sparse_matrix<representation, sparse::diagonal<2>, indexes, indexes>;

  const diagonal_matrix lhs{sparse::make_matrix<diagonal_matrix>(
      1. * m2, sparse::structural_zero, sparse::structural_zero, 2. * m2)};

  static_assert(std::same_as<decltype(lhs.at<0, 1>()), area>);
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
  assert((product.at<0, 1>() == 4. * m4));
  assert((product.at<1, 0>() == 10. * m4));
  assert((product.at<1, 1>() == 12. * m4));

  const auto squared_value{fcarouge::operator*(lhs, lhs)};

  static_assert(decltype(squared_value)::matrix::structure.count == 2);
  assert((squared_value.at<1, 1>() == 4. * m4));
  assert((squared_value.at<1, 0>() == 0. * m4));

  // A sum, a difference, a negation, a scaling, a transposition.
  const auto sum{lhs + rhs};

  assert((sum.at<0, 0>() == 4. * m2));
  assert((sum.at<0, 1>() == 4. * m2));
  assert(((lhs - rhs).at<1, 0>() == -5. * m2));
  assert(((-lhs).at<1, 1>() == -2. * m2));
  assert(((lhs * 2.).at<1, 1>() == 4. * m2));
  assert((transposed(rhs).at<0, 1>() == 5. * m2));

  // A norm in the element type.
  assert(matrix_one_norm(lhs) == 2. * m2);
  assert(matrix_inf_norm(rhs) == 11. * m2);

  // Heterogeneous elements with a structural zero of its own unit.
  using state =
      sparse_column_vector<representation, sparse::make_pattern<2, 1>({{0, 0}}),
                           length, duration_type>;

  state value{sparse::make_matrix<state>(2. * m, sparse::structural_zero)};

  static_assert(std::same_as<decltype(value.at<1>()), const duration_type>);
  assert(value.at<0>() == 2. * m);
  assert(value.at<1>() == 0. * s);

  value.at<0>(3. * m);

  assert(value.at<0>() == 3. * m);

  const auto [position, time_value]{value};

  assert(position == 3. * m);
  assert(time_value == 0. * s);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
