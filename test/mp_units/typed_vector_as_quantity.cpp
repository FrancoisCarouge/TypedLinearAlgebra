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
#include <format>

QUANTITY_SPEC(flight_velocity, mp_units::isq::velocity);
QUANTITY_SPEC(forward_velocity, mp_units::isq::velocity, mp_units::is_kind);
QUANTITY_SPEC(lateral_velocity, mp_units::isq::velocity, mp_units::is_kind);
QUANTITY_SPEC(vertical_velocity, mp_units::isq::velocity, mp_units::is_kind);

template <>
struct mp_units::vector_components<flight_velocity>
    : mp_units::vector_axes<forward_velocity, lateral_velocity,
                            vertical_velocity> {};

namespace fcarouge::test {
namespace {
// Inject the magnitude customization-point object (CPO, a variable with a call
// operator) so it competes with the library's magnitude free function.
using mp_units::magnitude;

using mp_units::si::unit_symbols::h;
using mp_units::si::unit_symbols::km;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;

//! @test Verifies typed linear algebra's mp-units compatibility plugin.
//!
//! @details Explore areas of interests.
[[maybe_unused]] const auto test{[] -> int {
  // The typed linear algebra library's vector can be used as a representation
  // for a quantity:
  {
    // A 3-component double-typed column vector, backed by an Eigen vector.
    using vector3d = column_vector<double, double, double, double>;

    // The homogeneous, double-typed vector satisfies mp-units's representation
    // requirements to back up a vector quantity.
    mp_units::quantity velocity{
        mp_units::isq::velocity(vector3d{30, 40, 0} * km / h)};

    // From there the framework treats the vector exactly like any other
    // representation:
    // * Unit conversion scales every component:
    mp_units::quantity velocity_mps{velocity.in(m / s)};
    // * vector quantity × scalar quantity:
    mp_units::quantity displacement{velocity * (2 * h)};
    // * Euclidean magnitude → scalar quantity:
    mp_units::quantity speed{magnitude(velocity)};

    assert(std::format("{}", velocity) == "[[30], [40], [0]] km/h");
    assert(std::format("{}", velocity_mps) ==
           "[[8.333333333333334], [11.11111111111111], [0]] m/s");
    assert(std::format("{}", displacement) == "[[60], [80], [0]] km");
    assert(std::format("{}", speed) == "50 km/h");
    assert(std::format("{}", velocity.numerical_value_in(velocity.unit)) ==
           "[[30], [40], [0]]");
  }

  // An mp-units vector quantity with strongly typed axes can be converted to a
  // heterogeneous typed vector of quantities.
  {
    // A 3-component double-typed column vector, backed by an Eigen vector.
    using vector3d = column_vector<double, double, double, double>;

    // A velocity vector quantity with strongly typed axes backed by a
    // 3-component double-typed column vector, backed by an Eigen vector.
    mp_units::quantity velocity{flight_velocity(vector3d{30, 40, 0} * km / h)};

    // A heterogeneous 3-component quantity-typed column vector, backed by an
    // Eigen vector.
    using velocities =
        column_vector<double, mp_units::quantity<forward_velocity[km / h]>,
                      mp_units::quantity<lateral_velocity[km / h]>,
                      mp_units::quantity<vertical_velocity[km / h]>>;

    // Convert the vector quantity to a heterogeneous typed column vector of
    // quantities, out of the box.
    velocities state{velocity};

    // Note the difference in presentation of the vector state compared to the
    // quantity vector:
    assert(std::format("{}", state) == "[[30 km/h], [40 km/h], [0 km/h]]");
    assert(std::format("{}", velocity) == "[[30], [40], [0]] km/h");

    // The vector quantity has a magnitude, but not the heterogeneous typed
    // vector of quantities, which doesn't have a meaning of a mathematical
    // vector nor homogeneous units.
    mp_units::quantity speed{magnitude(velocity)}; // 🗸
    // quantity speed{magnitude(state)}; // 𐄂

    assert(std::format("{}", speed) == "50 km/h");
  }

  // An mp-units vector quantity with strongly typed axes can be converted to a
  // homogeneous typed vector of quantities given compatible quantity
  // specifications.
  {
    // A 3-component double-typed column vector, backed by an Eigen vector.
    using vector3d = column_vector<double, double, double, double>;

    // A velocity vector quantity with strongly typed axes backed by a
    // 3-component double-typed column vector, backed by an Eigen vector.
    mp_units::quantity velocity{flight_velocity(vector3d{30, 40, 0} * km / h)};

    assert(std::format("{}", get<0>(velocity)) == "30 km/h");
    assert(std::format("{}", get<forward_velocity>(velocity)) == "30 km/h");

    // The library also permits the conversion to a strongly typed column
    // vector, backed by Eigen, with three components, each of type ISQ velocity
    // quantity, if the conversion supports the quantity specific cast.
    using velocities =
        column_vector<double,
                      mp_units::quantity<mp_units::isq::velocity[km / h]>,
                      mp_units::quantity<mp_units::isq::velocity[km / h]>,
                      mp_units::quantity<mp_units::isq::velocity[km / h]>>;

    // Convert the vector quantity to a homogeneous typed column vector of
    // quantities, through the plug-in tuple support.
    velocities state{velocity};

    // Dragons: magnitude possible in C++; that's not a mathematical vector.
    assert(std::format("{}", magnitude(state)) == "50 km/h");
  }

  // An mp-units vector quantity without strongly typed axes cannot be converted
  // to a homogeneous typed vector of quantities.
  {
    // A 3-component double-typed column vector, backed by an Eigen vector.
    using vector3d = column_vector<double, double, double, double>;

    // The homogeneous, double-typed vector satisfies mp-units's representation
    // requirements to back up a vector quantity.
    mp_units::quantity velocity{
        mp_units::isq::velocity(vector3d{30, 40, 0} * km / h)};

    assert(std::format("{}", velocity) == "[[30], [40], [0]] km/h");

    // Not supported:
    // assert(std::format("{}", get<0>(velocity)) == "30 km/h");

    using velocities =
        column_vector<double,
                      mp_units::quantity<mp_units::isq::velocity[km / h]>,
                      mp_units::quantity<mp_units::isq::velocity[km / h]>,
                      mp_units::quantity<mp_units::isq::velocity[km / h]>>;

    velocities state{velocity};

    assert(std::format("{}", state) == "[[30 km/h], [40 km/h], [0 km/h]]");
  }

  {
    // A 3-component double-typed column vector, backed by an Eigen vector.
    using vector3d = column_vector<double, double, double, double>;

    // The homogeneous, double-typed vector satisfies mp-units's representation
    // requirements to back up a vector quantity.
    mp_units::quantity velocity{flight_velocity(vector3d{30, 40, 0} * km / h)};

    [[maybe_unused]] auto [x, y, z] = velocity;
  }

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
