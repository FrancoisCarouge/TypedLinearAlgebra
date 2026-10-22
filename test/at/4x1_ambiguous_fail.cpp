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

namespace fcarouge::test {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;

using position = quantity<mp_units::isq::length[m]>;
using duration = quantity<mp_units::isq::duration[s]>;

//! @test The by-type `at` accessor rejects, at compile time, a request whose
//! type is the element type of several positions of a vector, even though
//! another, unrepeated element type remains locatable.
[[maybe_unused]] const auto test{[] -> int {
  column_vector<representation, position, position, position, duration> x{
      1. * m, 2. * m, 3. * m, 4. * s};

  // Intended: [[maybe_unused]] const auto value{x.at<duration>()};
  [[maybe_unused]] const auto value{x.at<position>()};

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
