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

#include <cassert>
#include <chrono>
#include <cstddef>
#include <mdspan>

namespace fcarouge::test {
using literals::operator""_i;
using representation = double;

namespace {
using seconds = std::chrono::duration<representation>;
using instant = std::chrono::time_point<std::chrono::steady_clock, seconds>;

//! @test The by-type `at` accessor over a distinct row vector backed by
//! std::mdspan storage. A duration and a time point neither convert to one
//! another nor share a common type, so the pair stays distinct.
[[maybe_unused]] const auto test{[] -> int {
  representation storage[]{0., 0.};
  std::mdspan span{&storage[0], std::extents<std::size_t, 1, 2>{}};

  row_vector<representation, seconds, instant> x{span};

  x.at<0_i>(seconds{2.});
  x.at<1_i>(instant{seconds{3.}});

  assert(seconds{2.} == x.at<seconds>());
  assert(seconds{2.} == x.at<0>());
  assert(seconds{2.} == x.at<0_i>());

  assert(instant{seconds{3.}} == x.at<instant>());
  assert(instant{seconds{3.}} == x.at<1>());
  assert(instant{seconds{3.}} == x.at<1_i>());

  x.at<0_i>(seconds{5.});

  assert(seconds{5.} == x.at<seconds>());

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
