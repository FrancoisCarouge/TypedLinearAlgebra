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

#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <mdspan>
#include <ratio>

namespace fcarouge::test {
namespace {
//! @test Verifies the `add` function over row vectors of distinct std::chrono
//! duration periods viewing caller-owned `std::array` memory.
[[maybe_unused]] const auto test{[] -> int {
  using seconds = std::chrono::duration<double>;
  using minutes = std::chrono::duration<double, std::ratio<60>>;

  std::array storage_a{1., 2.};
  std::array storage_b{3., 4.};
  std::array<double, 2> storage_r{};

  std::mdspan span_a{storage_a.data(), std::extents<std::size_t, 1, 2>{}};
  std::mdspan span_b{storage_b.data(), std::extents<std::size_t, 1, 2>{}};
  std::mdspan span_r{storage_r.data(), std::extents<std::size_t, 1, 2>{}};

  row_vector<double, seconds, minutes> a{span_a};
  row_vector<double, seconds, minutes> b{span_b};
  row_vector<double, seconds, minutes> r{span_r};

  add(a, b, r);

  assert(a.at<0>() == seconds{1.});
  assert(a.at<1>() == minutes{2.});
  assert(r.at<0>() == seconds{4.});
  assert(r.at<1>() == minutes{6.});
  assert((storage_r == std::array{4., 6.}));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
