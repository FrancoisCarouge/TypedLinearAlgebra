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
//! @brief Unit safe sparse linear algebra with mp-units.
//!
//! @details Demonstrate the compile-time sparsity pattern backend: only the
//! stored elements occupy memory and cost arithmetic operations, while every
//! element remains strongly typed, a structural zero included. Every operation
//! is real-time safe, without memory allocation, and available at compile time.

#include "fcarouge/linalg.hpp"

#include <cassert>
#include <tuple>

namespace fcarouge::sample {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::one;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;
using literals::operator""_o;

using position = quantity<mp_units::isq::length[m]>;
using velocity = quantity<mp_units::isq::speed[m / s]>;
using inverse_position = quantity<one / mp_units::isq::length[m]>;
using inverse_velocity = quantity<one / mp_units::isq::speed[m / s]>;

//! @brief Strongly typed sparse linear algebra samples.
[[maybe_unused]] const auto sample{[] -> int {
  // A state transition, a linear map: the column indexes are reciprocal to the
  // state's. The pattern is part of the type: the bottom left element is a
  // structural zero, without storage.
  using transition =
      typed_sparse_matrix<representation, sparse::upper_triangular<2>,
                          std::tuple<position, velocity>,
                          std::tuple<inverse_position, inverse_velocity>>;

  static_assert(sizeof(transition) == 3 * sizeof(representation));

  // Each stored element is written with its own type: a ratio, a duration.
  transition model;

  model.at<0, 0>(1. * one);
  model.at<0, 1>(2. * s);
  model.at<1, 1>(1. * one);

  // A structural zero is read with its own type, an inverse duration, by value.
  // Writing it does not compile: `model.at<1, 0>(1. / s);`.
  assert((model.at<1, 0>() == 0. / s));

  // The product only computes the terms the patterns admit: three of the four
  // of a dense product of these shapes.
  using state = typed_sparse_column_vector<representation, sparse::dense<2, 1>,
                                           position, velocity>;

  const state initial{1. * m, 3. * m / s};
  const state predicted{model * initial};

  assert(predicted.at<0>() == 7. * m);
  assert(predicted.at<1>() == 3. * m / s);

  // A pattern of any shape is a list of positions.
  constexpr sparse::pattern ends{sparse::shape<3, 1>, {{0, 0}, {2, 0}}};

  using sparse_state = typed_sparse_column_vector<representation, ends,
                                                  position, velocity, position>;

  static_assert(sizeof(sparse_state) == 2 * sizeof(representation));

  // Each operation computes the pattern of its result at compile time.
  const sparse_state first{
      sparse::make_matrix<sparse_state>(1. * m, 0_o, 2. * m)};
  const auto sum{first + first};

  static_assert(decltype(sum)::matrix::structure == ends);
  assert(sum.at<2>() == 4. * m);

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
