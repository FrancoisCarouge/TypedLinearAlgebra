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

#include <mp-units/systems/isq.h>
#include <mp-units/systems/si.h>

#include <functional>
#include <tuple>

namespace fcarouge::test {
namespace {
using mp_units::si::unit_symbols::km;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using position = quantity<mp_units::isq::length[m]>;
using kilo_position = quantity<mp_units::isq::length[km]>;
using velocity = quantity<mp_units::isq::velocity[m / s]>;
using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

//! @test Verifies the distinct typed matrix concept.
[[maybe_unused]] const auto test{[] -> int {
  {
    //! @todo Permit non-tuple declaration.
    using sq = matrix<double, std::tuple<position>, std::tuple<std::identity>>;
    static_assert(distinct_typed_matrix<sq>);
    static_assert(uniform_typed_matrix<sq>); // that's funny
  }
  {
    using cv = matrix<double, std::tuple<position, velocity>,
                      std::tuple<std::identity>>;
    static_assert(distinct_typed_matrix<cv>);
    static_assert(not uniform_typed_matrix<cv>);
  }
  {
    using cv = matrix<double, std::tuple<velocity, velocity>,
                      std::tuple<std::identity>>;
    static_assert(not distinct_typed_matrix<cv>);
    static_assert(uniform_typed_matrix<cv>);
  }
  {
    using rv = matrix<double, std::tuple<std::identity>,
                      std::tuple<position, velocity>>;
    static_assert(distinct_typed_matrix<rv>);
    static_assert(not uniform_typed_matrix<rv>);
  }
  {
    using rv = matrix<double, std::tuple<std::identity>,
                      std::tuple<position, position>>;
    static_assert(not distinct_typed_matrix<rv>);
    static_assert(uniform_typed_matrix<rv>);
  }
  {
    using mx = matrix<double, std::tuple<double, position>,
                      std::tuple<velocity, acceleration>>;
    static_assert(distinct_typed_matrix<mx>);
    static_assert(not uniform_typed_matrix<mx>);
  }
  {
    using mx = matrix<double, std::tuple<position, velocity>,
                      std::tuple<velocity, acceleration>>;
    static_assert(not distinct_typed_matrix<mx>);
    static_assert(not uniform_typed_matrix<mx>);
  }
  {
    using mx = matrix<double, std::tuple<position, position>,
                      std::tuple<position, position>>;
    static_assert(not distinct_typed_matrix<mx>);
    static_assert(uniform_typed_matrix<mx>);
  }
  {
    // Convertible but not identical element types: neither distinct nor
    // uniform. This is the case that separates distinctness from mere
    // non-uniformity.
    using cv = matrix<double, std::tuple<position, kilo_position>,
                      std::tuple<std::identity>>;
    static_assert(not distinct_typed_matrix<cv>);
    static_assert(not uniform_typed_matrix<cv>);
  }
  {
    // Row vector wider than two: every pair of the three element types must
    // be checked.
    using rv = matrix<double, std::tuple<std::identity>,
                      std::tuple<position, velocity, acceleration>>;
    static_assert(distinct_typed_matrix<rv>);
    static_assert(not uniform_typed_matrix<rv>);
  }
  {
    // Row vector wider than two with one convertible pair among three.
    using rv = matrix<double, std::tuple<std::identity>,
                      std::tuple<position, velocity, kilo_position>>;
    static_assert(not distinct_typed_matrix<rv>);
    static_assert(not uniform_typed_matrix<rv>);
  }
  {
    // Non-square two-dimension matrix, all six element types distinct.
    using mx = matrix<double, std::tuple<double, position>,
                      std::tuple<velocity, acceleration, position>>;
    static_assert(distinct_typed_matrix<mx>);
    static_assert(not uniform_typed_matrix<mx>);
  }
  return 0;
}()};
} // namespace
} // namespace fcarouge::test
