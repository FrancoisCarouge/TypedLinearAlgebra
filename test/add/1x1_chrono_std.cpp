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
#include <cstddef>
#include <mdspan>

namespace fcarouge::test {
namespace {
//! @test Verifies the singleton by singleton matrix `add` function with
//! std::chrono element types.
[[maybe_unused]] const auto test{[] -> int {
  using seconds = std::chrono::duration<double>;

  double storage_a{0.};
  double storage_b{0.};
  double storage_r{0.};

  std::mdspan span_a{&storage_a, std::extents<std::size_t, 1, 1>{}};
  std::mdspan span_b{&storage_b, std::extents<std::size_t, 1, 1>{}};
  std::mdspan span_r{&storage_r, std::extents<std::size_t, 1, 1>{}};

  row_vector<double, seconds> a{span_a};
  row_vector<double, seconds> b{span_b};
  row_vector<double, seconds> r{span_r};

  a.at(seconds{2.});
  b.at(seconds{3.});
  add(a, b, r);

  assert(r.at() == seconds{5.});

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
