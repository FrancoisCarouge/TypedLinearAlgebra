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
#include <ratio>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the addition rejects vectors of different scales, seconds and
//! minutes: their stored numbers cannot be computed on directly.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = std::chrono::duration<representation>;
  using large_unit = std::chrono::duration<representation, std::ratio<60>>;
  using small_column = column_vector<representation, small_unit, small_unit>;
  using large_column = column_vector<representation, large_unit, large_unit>;

  const small_unit first{30.};
  const large_unit second{3.};
  const small_column lhs{first, first};
  const large_column rhs{second, second};

  // Intended: [[maybe_unused]] const auto result{lhs + small_column{rhs}};
  [[maybe_unused]] const auto result{lhs + rhs};

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
