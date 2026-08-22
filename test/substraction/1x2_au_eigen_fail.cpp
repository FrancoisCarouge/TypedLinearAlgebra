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

#include <au/units/meters.hh>
#include <au/units/seconds.hh>

#include <cassert>

namespace fcarouge::test {
using literals::operator""_i;
using representation = double;

namespace {
//! @test Verifies the substraction operator with non-trivial types.
[[maybe_unused]] const auto test{[] {
  using position = au::QuantityD<au::Meters>;
  using velocity = au::QuantityD<au::UnitQuotientT<au::Meters, au::Seconds>>;

  row_vector<representation, position, velocity> a{
      au::meters(1.), au::meters(2.) / au::seconds(1.)};

  // Intended:
  // row_vector<representation, position, velocity> b{au::meters(3.),
  //                                                   au::meters(1.) /
  //                                                       au::seconds(1.)};
  row_vector<representation, velocity, position> b{
      au::meters(3.) / au::seconds(1.), au::meters(1.)};

  row_vector<representation, position, velocity> r{a - b};

  assert(au::meters(-2.) == r.at<0_i>());
  assert(au::meters(-2.) == r.at<0>());
  assert(au::meters(-2.) == r[0_i]);
  assert(au::meters(-2.) == r(0_i));

  assert(au::meters(1.) / au::seconds(1.) == r.at<1_i>());
  assert(au::meters(1.) / au::seconds(1.) == r.at<1>());
  assert(au::meters(1.) / au::seconds(1.) == r[1_i]);
  assert(au::meters(1.) / au::seconds(1.) == r(1_i));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
