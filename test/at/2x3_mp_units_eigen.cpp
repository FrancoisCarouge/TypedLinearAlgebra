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

#include "fcarouge/linalg.hpp"

#include <cassert>
#include <functional>
#include <tuple>

#include <Eigen/Eigen>

namespace fcarouge::test {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::kg;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;

using length = quantity<mp_units::isq::length[m]>;
using mass = quantity<mp_units::isq::mass[kg]>;
using time = quantity<mp_units::isq::time[s]>;

//! @test Verifies the at member accessor looking up the element by type on a
//! distinct, non-square two-dimension matrix. A rectangular shape catches a
//! row/column transposition in the linear-index arithmetic that a square shape
//! cannot.
[[maybe_unused]] const auto test{[] -> int {
  using rows = std::tuple<length, mass>;
  using columns = std::tuple<std::identity, time, decltype(1. * s * s)>;
  using matrix =
      typed_matrix<Eigen::Matrix<representation, 2, 3>, rows, columns>;

  Eigen::Matrix<representation, 2, 3> storage;
  storage << 1., 2., 3., 4., 5., 6.;
  const matrix x{storage};

  // Row-major element types:
  //   (0,0) length     (0,1) length*time     (0,2) length*time^2
  //   (1,0) mass       (1,1) mass*time       (1,2) mass*time^2
  assert((x.at<length>() == x.at<0, 0>()));
  assert((x.at<decltype(1. * m * s)>() == x.at<0, 1>()));
  assert((x.at<decltype(1. * m * s * s)>() == x.at<0, 2>()));
  assert((x.at<mass>() == x.at<1, 0>()));
  assert((x.at<decltype(1. * kg * s)>() == x.at<1, 1>()));
  assert((x.at<decltype(1. * kg * s * s)>() == x.at<1, 2>()));

  assert((x.at<length>() == 1. * m));
  assert((x.at<decltype(1. * kg * s * s)>() == 6. * kg * s * s));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
