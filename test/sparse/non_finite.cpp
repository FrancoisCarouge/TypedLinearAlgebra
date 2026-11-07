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
#include <limits>

namespace fcarouge::test {
namespace {
//! @test Verifies the documented semantics of non-finite values: a structural
//! zero has no term, so it stays an exact zero, where a dense computation
//! yields a not-a-number value.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::diagonal;

  constexpr double infinity{std::numeric_limits<double>::infinity()};

  const sparse_matrix<double, diagonal<2>> value{
      sparse::make_matrix<sparse_matrix<double, diagonal<2>>>(
          1., sparse::structural_zero, sparse::structural_zero, 2.)};

  // Scaling.
  const auto scaled{value * infinity};

  assert((scaled.at<0, 0>() == infinity));
  assert((scaled.at<0, 1>() == 0.));

  const auto dense_scaled{sparse_matrix<double, sparse::dense<2, 2>>{value} *
                          infinity};

  assert((dense_scaled.at<0, 1>() != dense_scaled.at<0, 1>()));

  // Dividing by zero.
  const auto divided{value / 0.};

  assert((divided.at<0, 0>() == infinity));
  assert((divided.at<0, 1>() == 0.));

  // A product.
  const sparse_matrix<double, diagonal<2>> unbounded{
      sparse::make_matrix<sparse_matrix<double, diagonal<2>>>(
          infinity, sparse::structural_zero, sparse::structural_zero, 1.)};
  const auto product{unbounded * value};

  assert((product.at<0, 0>() == infinity));
  assert((product.at<0, 1>() == 0.));
  assert((product.at<1, 0>() == 0.));

  // A stored not-a-number value stays visible.
  constexpr double not_a_number{std::numeric_limits<double>::quiet_NaN()};

  const sparse_matrix<double, diagonal<2>> invalid{
      sparse::make_matrix<sparse_matrix<double, diagonal<2>>>(
          not_a_number, sparse::structural_zero, sparse::structural_zero, 1.)};

  assert((invalid != invalid));
  assert((matrix_frob_norm(invalid) != matrix_frob_norm(invalid)));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
