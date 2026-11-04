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
of the public at large_unit and to the detriment of our heirs and
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

#include <chrono>
#include <cstddef>
#include <mdspan>
#include <ratio>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the dot product rejects terms of different scales,
//! seconds and minutes.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = std::chrono::duration<representation>;
  using large_unit = std::chrono::duration<representation, std::ratio<60>>;

  representation storage_values[2]{};
  representation storage_ones[2]{};

  std::mdspan span_values{&storage_values[0],
                          std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_ones{&storage_ones[0], std::extents<std::size_t, 2, 1>{}};

  column_vector<representation, small_unit, large_unit> values{span_values};
  column_vector<representation, representation, representation> ones{span_ones};

  // Intended, on a uniform vector: dot(uniform_values, ones);
  [[maybe_unused]] const auto result{dot(values, ones)};

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
