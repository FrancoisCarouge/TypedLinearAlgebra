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
#include <concepts>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies adding a duration vector and a time point vector deduces a
//! time point result vector regardless of operand order, since std::chrono
//! defines both time point + duration and duration + time point as a time
//! point.
[[maybe_unused]] const auto test{[] -> int {
  using seconds = std::chrono::duration<representation>;
  using instant = std::chrono::time_point<std::chrono::system_clock, seconds>;

  const row_vector<representation, seconds, seconds> a{seconds{5.},
                                                       seconds{9.}};
  const row_vector<representation, instant, instant> b{instant{seconds{2.}},
                                                       instant{seconds{3.}}};
  const auto r{a + b};

  assert(instant{seconds{7.}} == r.at<0>());
  assert(instant{seconds{12.}} == r.at<1>());
  static_assert(
      std::same_as<decltype(r)::column_indexes, std::tuple<instant, instant>>);

  const auto s{b + a};

  assert(instant{seconds{7.}} == s.at<0>());
  assert(instant{seconds{12.}} == s.at<1>());
  static_assert(
      std::same_as<decltype(s)::column_indexes, std::tuple<instant, instant>>);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
