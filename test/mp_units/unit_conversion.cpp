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

#include <array>
#include <cassert>
#include <format>
#include <tuple>
#include <type_traits>

inline constexpr struct flight_velocity
    : mp_units::quantity_spec<mp_units::isq::velocity> {
} flight_velocity;
inline constexpr struct forward_velocity
    : mp_units::quantity_spec<mp_units::isq::velocity, mp_units::is_kind> {
} forward_velocity;
inline constexpr struct lateral_velocity
    : mp_units::quantity_spec<mp_units::isq::velocity, mp_units::is_kind> {
} lateral_velocity;
inline constexpr struct vertical_velocity
    : mp_units::quantity_spec<mp_units::isq::velocity, mp_units::is_kind> {
} vertical_velocity;

template <>
struct mp_units::component_vectors<flight_velocity>
    : mp_units::vector_axes<forward_velocity, lateral_velocity,
                            vertical_velocity> {};

namespace fcarouge::test {
namespace {
using mp_units::isq::length;
using mp_units::isq::speed;
using mp_units::isq::velocity;
using mp_units::si::unit_symbols::h;
using mp_units::si::unit_symbols::km;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;

// A 3-component double-typed column vector.
using vector3d = column_vector<double, double, double, double>;

// The per-axis flight velocity elements, declared in kilometers per hour.
using forward_kmph = mp_units::quantity<forward_velocity[km / h]>;
using lateral_kmph = mp_units::quantity<lateral_velocity[km / h]>;
using vertical_kmph = mp_units::quantity<vertical_velocity[km / h]>;
using velocities =
    column_vector<double, forward_kmph, lateral_kmph, vertical_kmph>;

// The same axes, swapped.
using swapped_velocities =
    column_vector<double, lateral_kmph, forward_kmph, vertical_kmph>;

// The per-axis flight velocity sources, in meters per second.
using forward_mps = mp_units::quantity<forward_velocity[m / s]>;
using lateral_mps = mp_units::quantity<lateral_velocity[m / s]>;
using vertical_mps = mp_units::quantity<vertical_velocity[m / s]>;

// A uniformly typed speed vector, declared in kilometers per hour.
using speed_kmph = mp_units::quantity<speed[km / h]>;
using speed_mps = mp_units::quantity<speed[m / s]>;
using speeds = column_vector<double, speed_kmph, speed_kmph>;

// The flight velocity components without their axes, declared in kilometers
// per hour.
using velocity_kmph = mp_units::quantity<velocity[km / h]>;
using velocity_mps = mp_units::quantity<velocity[m / s]>;
using generic_velocities =
    column_vector<double, velocity_kmph, velocity_kmph, velocity_kmph>;

// The flight velocity vector quantities.
using flight_velocity_kmph =
    mp_units::quantity<flight_velocity[km / h], vector3d>;
using flight_velocity_mps =
    mp_units::quantity<flight_velocity[m / s], vector3d>;

// A tuple-like value implicitly converts when each of its elements implicitly
// converts to the corresponding element type: same axes, any unit.
static_assert(std::is_constructible_v<
              velocities, std::tuple<forward_mps, lateral_mps, vertical_mps>>);
static_assert(std::is_constructible_v<velocities, flight_velocity_mps>);
static_assert(std::convertible_to<flight_velocity_mps, velocities>);
static_assert(other_tuple_like_vector<
              std::tuple<forward_mps, lateral_mps, vertical_mps>, velocities>);

// The axes are not interchangeable, neither from a tuple nor from a vector
// quantity decomposing in a different axes order.
static_assert(
    not std::is_constructible_v<
        velocities, std::tuple<lateral_kmph, forward_kmph, vertical_kmph>>);
static_assert(
    not std::is_constructible_v<swapped_velocities, flight_velocity_kmph>);
static_assert(
    not other_tuple_like_vector<
        std::tuple<lateral_kmph, forward_kmph, vertical_kmph>, velocities>);

// An axis and its more general quantity only explicitly convert to each other,
// so does a tuple-like value of them.
static_assert(std::is_constructible_v<forward_kmph, velocity_mps>);
static_assert(not std::is_convertible_v<velocity_mps, forward_kmph>);
static_assert(
    not tla::is_convertible_tuple_like<
        std::tuple<velocity_mps, velocity_mps, velocity_mps>, velocities>());
static_assert(
    std::is_constructible_v<
        velocities, std::tuple<velocity_mps, velocity_mps, velocity_mps>>);
static_assert(
    not std::convertible_to<
        std::tuple<velocity_mps, velocity_mps, velocity_mps>, velocities>);
static_assert(std::is_constructible_v<generic_velocities, flight_velocity_mps>);
static_assert(not std::convertible_to<flight_velocity_mps, generic_velocities>);

// Neither do unrelated quantities nor plain numbers.
using length_m = mp_units::quantity<length[m]>;
static_assert(not std::is_constructible_v<
              velocities, std::tuple<length_m, length_m, length_m>>);
static_assert(not std::is_constructible_v<velocities,
                                          std::tuple<double, double, double>>);
static_assert(not std::is_constructible_v<velocities, std::array<double, 3>>);

// Elements convert to a typed vector as a tuple-like value of them would: of
// the vector's size, implicitly when every element implicitly converts,
// explicitly when some element only explicitly converts.
static_assert(std::is_constructible_v<velocities, forward_mps, lateral_mps,
                                      vertical_mps>);
static_assert(
    not std::is_constructible_v<velocities, forward_mps, lateral_mps>);
static_assert(not std::is_constructible_v<velocities, forward_mps, lateral_mps,
                                          vertical_mps, vertical_mps>);
static_assert(not std::is_constructible_v<velocities, lateral_kmph,
                                          forward_kmph, vertical_kmph>);
static_assert(std::is_constructible_v<velocities, velocity_mps, velocity_mps,
                                      velocity_mps>);
static_assert(
    not std::is_constructible_v<velocities, length_m, length_m, length_m>);
static_assert(not std::is_constructible_v<velocities, double, double, double>);

//! @test Verifies the elements are converted to their declared unit, whatever
//! the construction path.
[[maybe_unused]] const auto test{[] -> int {
  const not_realtime opt_out;

  // The element caster stores the numerical value in the element's unit.
  {
    assert((cast<double, speed_kmph>(10. * speed[m / s]) == 36.));
  }

  // From a vector quantity in another unit.
  {
    const flight_velocity_mps value{
        flight_velocity(vector3d{10, -20, 5} * m / s)};
    const velocities state = value;

    assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  }

  // From a vector quantity in another unit, explicitly dropping its axes.
  {
    const flight_velocity_mps value{
        flight_velocity(vector3d{10, -20, 5} * m / s)};
    const generic_velocities state{value};

    assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  }

  // From a tuple of quantities in another unit.
  {
    const std::tuple value{forward_velocity(10. * m / s),
                           lateral_velocity(-20. * m / s),
                           vertical_velocity(5. * m / s)};
    const velocities state = value;

    assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  }

  // From a tuple of more general quantities in another unit, explicitly
  // specified to the axes.
  {
    const std::tuple value{velocity(10. * m / s), velocity(-20. * m / s),
                           velocity(5. * m / s)};
    const velocities state{value};

    assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  }

  // From an array of quantities in another unit.
  {
    const std::array<speed_mps, 2> values{10. * speed[m / s],
                                          20. * speed[m / s]};
    const speeds state = values;

    assert(std::format("{}", state) == "[[36 km/h], [72 km/h]]");
  }

  // From elements in another unit.
  {
    const velocities state{forward_velocity(10. * m / s),
                           lateral_velocity(-20. * m / s),
                           vertical_velocity(5. * m / s)};

    assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  }

  // From implicitly converted elements in another unit.
  {
    const velocities state = {forward_velocity(10. * m / s),
                              lateral_velocity(-20. * m / s),
                              vertical_velocity(5. * m / s)};

    assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  }

  // From more general elements in another unit, explicitly specified to the
  // axes.
  {
    const velocities state{velocity(10. * m / s), velocity(-20. * m / s),
                           velocity(5. * m / s)};

    assert(std::format("{}", state) == "[[36 km/h], [-72 km/h], [18 km/h]]");
  }

  // From list-initializers in another unit.
  {
    const speeds state{{10. * speed[m / s]}, {20. * speed[m / s]}};

    assert(std::format("{}", state) == "[[36 km/h], [72 km/h]]");
  }

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
