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
#include "fcarouge/realtime.hpp"

#include <cassert>
#include <format>

namespace fcarouge::test {
namespace {
//! @test Verifies a structural zero is formatted as an element, so the format
//! does not depend on the pattern.
[[maybe_unused]] const auto test{[] -> int {
  const not_realtime opt_out;

  using sparse::make_pattern;

  const sparse_matrix<double, sparse::diagonal<2>> value{
      sparse::make_matrix<sparse_matrix<double, sparse::diagonal<2>>>(
          1., sparse::structural_zero, sparse::structural_zero, 2.5)};
  const sparse_matrix<double, sparse::dense<2, 2>> dense_value{{1., 0.},
                                                               {0., 2.5}};

  assert(std::format("{}", value) == "[[1, 0], [0, 2.5]]");
  assert(std::format("{}", value) == std::format("{}", dense_value));

  const sparse_matrix<double, make_pattern<3, 1>({{1, 0}})> column{
      sparse::make_matrix<sparse_matrix<double, make_pattern<3, 1>({{1, 0}})>>(
          sparse::structural_zero, 4., sparse::structural_zero)};
  const sparse_matrix<double, make_pattern<1, 3>({{0, 2}})> row{
      sparse::make_matrix<sparse_matrix<double, make_pattern<1, 3>({{0, 2}})>>(
          sparse::structural_zero, sparse::structural_zero, 5.)};

  assert(std::format("{}", column) == "[[0], [4], [0]]");
  assert(std::format("{}", row) == "[0, 0, 5]");

  const sparse_matrix<double, sparse::zero<1, 1>> singleton;

  assert(std::format("{}", singleton) == "0");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
