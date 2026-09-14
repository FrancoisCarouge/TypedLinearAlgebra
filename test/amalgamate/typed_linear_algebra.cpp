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

//! @file
//! @brief Verifies the `amalgamate/fcarouge/typed_linear_algebra.h`
//! single-header amalgamation, generated for tools such as Compiler Explorer
//! that cannot resolve this project's multi-file layout, remains a
//! functional drop-in replacement for the multi-file headers. Eigen is used
//! here only as a readily available element storage backend with its own
//! built-in arithmetic; the amalgamation under test contains none of it.

#include "typed_linear_algebra.h"

#include <cassert>
#include <format>
#include <tuple>

#include <Eigen/Eigen>

namespace fcarouge::test {
namespace {
template <typename Representation, typename RowIndexes, typename ColumnIndexes>
using matrix =
    typed_matrix<Eigen::Matrix<Representation, std::tuple_size_v<RowIndexes>,
                               std::tuple_size_v<ColumnIndexes>>,
                 RowIndexes, ColumnIndexes>;

[[maybe_unused]] const auto test{[] -> int {
  using indexes = std::tuple<int, int, int>;
  const matrix<int, indexes, indexes> m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  assert(std::format("{}", m) == "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]");

  const matrix<int, indexes, indexes> n{m + m};
  assert(std::format("{}", n) == "[[2, 4, 6], [8, 10, 12], [14, 16, 18]]");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
