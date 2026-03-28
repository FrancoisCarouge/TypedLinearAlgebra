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

#include "fcarouge/linalg.hpp"

#include <cassert>

namespace fcarouge::test {
using literals::operator""_i;

namespace {
//! @test Verifies the integral literal operator accessor.
[[maybe_unused]] auto test{[] {
  matrix<double, 2, 2> m{{1., 2.}, {3., 4.}};

  assert((m(0_i, 0_i) == 1.));
  assert((m(0_i, 1_i) == 2.));
  assert((m(1_i, 0_i) == 3.));
  assert((m(1_i, 1_i) == 4.));

  m(0_i, 0_i) = 11.;
  m(0_i, 1_i) = 12.;
  m(1_i, 0_i) = 13.;
  m(1_i, 1_i) = 14.;

  assert((m(0_i, 0_i) == 11.));
  assert((m(0_i, 1_i) == 12.));
  assert((m(1_i, 0_i) == 13.));
  assert((m(1_i, 1_i) == 14.));

  m[0_i, 0_i] = 1.;
  m[0_i, 1_i] = 2.;
  m[1_i, 0_i] = 3.;
  m[1_i, 1_i] = 4.;

  assert((m[0_i, 0_i] == 1.));
  assert((m[0_i, 1_i] == 2.));
  assert((m[1_i, 0_i] == 3.));
  assert((m[1_i, 1_i] == 4.));

  m.at<0_i, 0_i>() = 11.;
  m.at<0_i, 1_i>() = 12.;
  m.at<1_i, 0_i>() = 13.;
  m.at<1_i, 1_i>() = 14.;

  assert((m.at<0_i, 0_i>() == 11.));
  assert((m.at<0_i, 1_i>() == 12.));
  assert((m.at<1_i, 0_i>() == 13.));
  assert((m.at<1_i, 1_i>() == 14.));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
