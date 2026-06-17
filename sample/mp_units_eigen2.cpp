/* Typed Linear Algebra
Version 0.2.0
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

namespace fcarouge::sample {
namespace {
template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, double>;

template <typename... Types>
using column_vector =
    typed_column_vector<Eigen::Vector<double, sizeof...(Types)>, Types...>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;
using velocity = quantity<mp_units::isq::velocity[m / s]>;

// Vector of homogeneous quantities [velocity, velocity, velocity].
using vector_of_3_velocities = column_vector<velocity, velocity, velocity>;

// 3-dimensional vector of velocity quantity of velocities.
using velocity_3d_vector =
    mp_units::quantity<mp_units::isq::velocity[m / s], vector_of_3_velocities>;

// Future: [V, V, V] is not as safe/nice as [Vx, Vy, Vz]!

[[maybe_unused]] const auto sample{[] {
  // Look up ^^

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
