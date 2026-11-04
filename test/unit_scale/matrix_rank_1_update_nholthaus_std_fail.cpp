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

#include <units/length.h>

#include <cstddef>
#include <mdspan>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the rank-1 update rejects a result of a different scale than
//! its products, meters and kilometers.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = units::length::meters<representation>;
  using large_unit = units::length::kilometers<representation>;
  using scalars = std::tuple<representation, representation>;

  representation storage_lhs[2]{};
  representation storage_rhs[2]{};
  representation storage_result[4]{};

  std::mdspan span_lhs{&storage_lhs[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_rhs{&storage_rhs[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_result{&storage_result[0],
                          std::extents<std::size_t, 2, 2>{}};

  column_vector<representation, small_unit, large_unit> lhs{span_lhs};
  column_vector<representation, representation, representation> rhs{span_rhs};
  matrix<representation, std::tuple<small_unit, small_unit>, scalars> result{
      span_result};

  // Intended, on a result of the same units: matrix_rank_1_update(lhs, rhs,
  // mixed);
  matrix_rank_1_update(lhs, rhs, result);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
