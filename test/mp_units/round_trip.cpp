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
#include "fcarouge/realtime.hpp"

#include <cassert>
#include <concepts>
#include <format>
#include <type_traits>

QUANTITY_SPEC(flight_velocity, mp_units::isq::velocity);
QUANTITY_SPEC(forward_velocity, mp_units::isq::velocity, mp_units::is_kind);
QUANTITY_SPEC(lateral_velocity, mp_units::isq::velocity, mp_units::is_kind);
QUANTITY_SPEC(vertical_velocity, mp_units::isq::velocity, mp_units::is_kind);

template <>
struct mp_units::component_vectors<flight_velocity>
    : mp_units::vector_axes<forward_velocity, lateral_velocity,
                            vertical_velocity> {};

namespace fcarouge::test {
namespace {
using mp_units::si::unit_symbols::h;
using mp_units::si::unit_symbols::km;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;

// A 3-component double-typed column vector, backed by an Eigen vector.
using vector3d = column_vector<double, double, double, double>;

// A heterogeneous 3-component quantity-typed column vector, one quantity per
// named axis of the flight velocity, backed by an Eigen vector.
using velocities =
    column_vector<double, mp_units::quantity<forward_velocity[km / h]>,
                  mp_units::quantity<lateral_velocity[km / h]>,
                  mp_units::quantity<vertical_velocity[km / h]>>;

// The same axes, swapped.
using swapped_velocities =
    column_vector<double, mp_units::quantity<lateral_velocity[km / h]>,
                  mp_units::quantity<forward_velocity[km / h]>,
                  mp_units::quantity<vertical_velocity[km / h]>>;

// The flight velocity vector quantities, backed by the double-typed vector.
using flight_velocity_kmph =
    mp_units::quantity<flight_velocity[km / h], vector3d>;
using flight_velocity_mps =
    mp_units::quantity<flight_velocity[m / s], vector3d>;

// A typed matrix converts to another only element by element: a quantity does
// not implicitly convert to a double.
static_assert(not std::convertible_to<velocities, vector3d>);
static_assert(not std::constructible_from<vector3d, velocities>);
static_assert(not std::is_assignable_v<vector3d &, velocities>);
static_assert(not std::is_assignable_v<vector3d &, const velocities &>);
static_assert(not convertible_typed_matrix<velocities, vector3d>);

// Element-wise compatible conversions remain advertised.
static_assert(std::convertible_to<velocities, velocities>);
static_assert(std::convertible_to<vector3d, vector3d>);
static_assert(convertible_typed_matrix<velocities, velocities>);

// A vector quantity composes from the heterogeneous typed vector, in the axes'
// order only.
static_assert(std::convertible_to<velocities, flight_velocity_kmph>);
static_assert(std::convertible_to<velocities, flight_velocity_mps>);
static_assert(
    not std::constructible_from<flight_velocity_kmph, swapped_velocities>);

//! @test Verifies a vector quantity with strongly typed axes round trips
//! through a heterogeneous typed vector of quantities.
[[maybe_unused]] const auto test{[] -> int {
  const not_realtime opt_out;
  const flight_velocity_kmph velocity{
      flight_velocity(vector3d{36, -72, 18} * km / h)};

  // Decompose the vector quantity into the heterogeneous typed vector.
  const velocities state{velocity};

  // Compose a vector quantity back from the heterogeneous typed vector, in the
  // same unit and in another unit.
  const flight_velocity_kmph round_trip{state};
  const flight_velocity_mps converted{state};

  assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  assert(round_trip == velocity);
  assert(std::format("{}", round_trip) == "[[36], [-72], [18]] km/h");
  assert(converted == velocity);
  assert(std::format("{}", converted) == "[[10], [-20], [5]] m/s");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
