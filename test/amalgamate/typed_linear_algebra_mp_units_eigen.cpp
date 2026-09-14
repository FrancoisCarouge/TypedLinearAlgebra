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
//! @brief Verifies the
//! `amalgamate/fcarouge/typed_linear_algebra_mp_units_eigen.h` single-header
//! amalgamation, generated for tools such as Compiler Explorer that cannot
//! resolve this project's multi-file layout, remains a functional drop-in
//! replacement for the multi-file `fcarouge/linalg.hpp` mp-units and Eigen
//! backend. An excerpt of the assertions from `sample/mp_units_eigen.cpp`.

#include "typed_linear_algebra_mp_units_eigen.h"

#include <cassert>
#include <format>

namespace fcarouge::test {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

using position = quantity<mp_units::isq::length[m]>;
using velocity = quantity<mp_units::isq::velocity[m / s]>;
using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

template <typename... Types>
using column_vector = column_vector<representation, Types...>;

[[maybe_unused]] const auto test{[] -> int {
  using state = column_vector<position, velocity, acceleration>;

  state x0{3. * m, 2. * m / s, 1. * m / s2};

  assert(std::format("{}", x0) == "[[3 m], [2 m/s], [1 m/s²]]");

  x0.at<1>(2.5 * m / s);
  assert(x0.at<1>() == 2.5 * m / s);

  state x1{x0 * 3.};
  assert(std::format("{}", x1) == "[[9 m], [7.5 m/s], [3 m/s²]]");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
