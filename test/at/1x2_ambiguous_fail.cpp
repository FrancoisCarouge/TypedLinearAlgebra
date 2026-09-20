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

#include <functional>
#include <tuple>

#include <Eigen/Eigen>

namespace fcarouge::test {
namespace {
// Two positions publicly inheriting a common, otherwise unrelated base: the
// distinctness check does not see the shared base as a conversion target, so
// the matrix stays distinct (see distinct_typed_matrix/mp_units_eigen.cpp's
// known over-approximation), yet both positions convert to it.
struct base {};
struct first : base {};
struct second : base {};

using row = std::tuple<std::identity>;
using column = std::tuple<first, second>;

//! @test The by-type `at` accessor rejects, at compile time, a request whose
//! type more than one element converts to, even on an otherwise distinct
//! matrix.
[[maybe_unused]] const auto test{[] -> int {
  typed_matrix<Eigen::Matrix<double, 1, 2>, row, column> x{};

  [[maybe_unused]] const auto value{x.at<base>()};

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
