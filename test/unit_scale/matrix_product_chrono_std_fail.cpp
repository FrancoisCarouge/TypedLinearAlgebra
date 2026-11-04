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
#include <functional>
#include <mdspan>
#include <ratio>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the matrix product function rejects a result of a different
//! scale than its terms, seconds and minutes.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = std::chrono::duration<representation>;
  using large_unit = std::chrono::duration<representation, std::ratio<60>>;
  using scalars = std::tuple<representation, representation>;
  using gain =
      matrix<representation, scalars, std::tuple<std::identity, std::identity>>;
  using mixed_rows =
      matrix<representation, std::tuple<small_unit, large_unit>, scalars>;
  using small_rows =
      matrix<representation, std::tuple<small_unit, small_unit>, scalars>;

  representation storage_lhs[4]{};
  representation storage_identity[4]{};
  representation storage_result[4]{};

  std::mdspan span_lhs{&storage_lhs[0], std::extents<std::size_t, 2, 2>{}};
  std::mdspan span_identity{&storage_identity[0],
                            std::extents<std::size_t, 2, 2>{}};
  std::mdspan span_result{&storage_result[0],
                          std::extents<std::size_t, 2, 2>{}};

  mixed_rows lhs{span_lhs};
  gain identity{span_identity};
  small_rows result{span_result};

  // Intended, on a result of the same units: matrix_product(lhs, identity,
  // mixed);
  matrix_product(lhs, identity, result);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
