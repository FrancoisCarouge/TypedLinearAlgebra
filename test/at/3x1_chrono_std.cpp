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

//! @test The compile-time and runtime indexed accessors over a column vector
//! backed by a column-major, `std::layout_left` std::mdspan storage. The
//! vector's linear index addresses the storage's row: reading it as the
//! storage's column instead lands, through the column-major mapping, on the
//! trailing values of the oversized buffer rather than on the vector's own.
[[maybe_unused]] const auto test{[] -> int {
  [[maybe_unused]] const std::size_t index{1};

  representation storage[]{1., 2., 3., 40., 50., 60., 70.};
  std::mdspan<representation, std::extents<std::size_t, 3, 1>, std::layout_left>
      span{&storage[0]};
  typed_column_vector<decltype(span), seconds, seconds, seconds> vector{span};

  assert(seconds{1.} == vector.at<0>());
  assert(seconds{2.} == vector.at<1_i>());
  assert(seconds{3.} == vector.at<2>());
  assert(seconds{2.} == vector[index]);
  assert(seconds{2.} == vector(index));

  vector.at<2_i>(seconds{9.});

  assert(seconds{9.} == vector.at<2>());
  assert(9. == storage[2]);
  assert(70. == storage[6]);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
