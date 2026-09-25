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

#include <au/units/meters.hh>
#include <au/units/seconds.hh>

#include <cassert>
#include <concepts>

namespace fcarouge::test {
namespace {
using representation = double;

//! @test The by-type `at` accessor resolves each element of a distinct row
//! vector to the same position, value, and type its integral index does, Au
//! quantities, Eigen-family backend.
[[maybe_unused]] const auto test{[] -> int {
  using au::symbols::m;
  using au::symbols::s;

  using position = au::QuantityD<au::Meters>;
  using velocity = au::QuantityD<au::UnitQuotientT<au::Meters, au::Seconds>>;

  row_vector<representation, position, velocity> x{2. * m, 3. * m / s};

  assert(x.at<position>() == 2. * m);
  assert(x.at<velocity>() == 3. * m / s);

  assert(x.at<position>() == x.at<0>());
  assert(x.at<velocity>() == x.at<1>());

  static_assert(std::same_as<decltype(x.at<velocity>()), decltype(x.at<1>())>);

  x.at<1>(9. * m / s);

  assert(x.at<velocity>() == 9. * m / s);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
