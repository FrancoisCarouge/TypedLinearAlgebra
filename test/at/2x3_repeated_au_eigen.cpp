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

#include <au/units/grams.hh>
#include <au/units/meters.hh>
#include <au/units/seconds.hh>

#include <cassert>
#include <concepts>
#include <functional>
#include <tuple>

namespace fcarouge::test {
namespace {
using representation = double;

template <typename Matrix, typename Type>
concept locatable = requires(Matrix value) { value.template at<Type>(); };

//! @test The by-type `at` accessor resolves the elements of a non-uniform,
//! non-square two-dimension matrix whose types occur once, while the types
//! repeated across several positions are not locatable, Au quantities,
//! Eigen-family backend.
[[maybe_unused]] const auto test{[] -> int {
  using au::symbols::g;
  using au::symbols::m;
  using au::symbols::s;

  using length = au::QuantityD<au::Meters>;
  using mass = au::QuantityD<au::Grams>;
  using time = au::QuantityD<au::Seconds>;
  using length_time = au::QuantityD<au::UnitProductT<au::Meters, au::Seconds>>;
  using mass_time = au::QuantityD<au::UnitProductT<au::Grams, au::Seconds>>;

  using rows = std::tuple<length, mass>;
  using columns = std::tuple<std::identity, std::identity, time>;
  using state = matrix<representation, rows, columns>;

  // Row-major element types:
  //   (0,0) length   (0,1) length   (0,2) length*time
  //   (1,0) mass     (1,1) mass     (1,2) mass*time
  state x;

  x.at<0, 0>(1. * m);
  x.at<0, 1>(2. * m);
  x.at<0, 2>(3. * m * s);
  x.at<1, 0>(4. * g);
  x.at<1, 1>(5. * g);
  x.at<1, 2>(6. * g * s);

  assert((x.at<length_time>() == 3. * m * s));
  assert((x.at<mass_time>() == 6. * g * s));
  assert((x.at<length_time>() == x.at<0, 2>()));
  assert((x.at<mass_time>() == x.at<1, 2>()));

  static_assert(
      std::same_as<decltype(x.at<length_time>()), decltype(x.at<0, 2>())>);
  static_assert(
      std::same_as<decltype(x.at<mass_time>()), decltype(x.at<1, 2>())>);
  static_assert(not uniform_typed_matrix<state>);
  static_assert(locatable<state, length_time>);
  static_assert(locatable<state, mass_time>);
  static_assert(not locatable<state, length>);
  static_assert(not locatable<state, mass>);
  static_assert(not locatable<state, time>);

  x.at<1, 2>(9. * g * s);

  assert((x.at<mass_time>() == 9. * g * s));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
