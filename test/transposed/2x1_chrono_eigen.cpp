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
#include <chrono>
#include <concepts>
#include <ratio>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the transposed algorithm turns a column vector of
//! std::chrono durations of distinct periods into the matching row vector,
//! Eigen backend.
[[maybe_unused]] const auto test{[] -> int {
  using seconds = std::chrono::duration<representation>;
  using minutes = std::chrono::duration<representation, std::ratio<60>>;

  column_vector<representation, seconds, minutes> n{seconds{42.}, minutes{43.}};

  assert(n.at<0>() == seconds{42.});
  assert(n.at<1>() == minutes{43.});
  static_assert(
      std::same_as<decltype(n)::row_indexes, std::tuple<seconds, minutes>>);

  row_vector<representation, seconds, minutes> nᵀ{transposed(n)};

  assert(nᵀ.at<0>() == seconds{42.});
  assert(nᵀ.at<1>() == minutes{43.});
  static_assert(
      std::same_as<decltype(nᵀ)::column_indexes, std::tuple<seconds, minutes>>);

  n.at<0>(seconds{24.});
  assert(nᵀ.at<0>() == seconds{42.});

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
