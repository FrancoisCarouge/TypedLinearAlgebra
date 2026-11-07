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

namespace fcarouge::test {
namespace {
//! @test Verifies every operation at run time, in the real-time context of the
//! real-time sanitizer job: no memory allocation, no lock, no blocking call.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::diagonal;
  using sparse::make_pattern;
  using sparse::tridiagonal;

  using band_matrix = sparse_matrix<double, tridiagonal<3>>;

  sparse_matrix<double, diagonal<3>> scale_matrix;

  scale_matrix.at<0, 0>(2.);
  scale_matrix.at<1, 1>(3.);
  scale_matrix.at<2, 2>(4.);

  const band_matrix band{
      sparse::make_matrix<band_matrix>(1., 2., sparse::structural_zero, 3., 4.,
                                       5., sparse::structural_zero, 6., 7.)};
  const band_matrix copy{band};

  const auto product{scale_matrix * band};
  const auto sum{scale_matrix + band};
  const auto difference{scale_matrix - band};
  const auto negated{-band};
  const auto scaled{band * 2.};
  const auto divided{band / 2.};
  const auto mirrored{transposed(band)};
  const sparse_matrix<double, sparse::dense<3, 3>> widened{band};

  const sparse_matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})> column{
      sparse::make_matrix<
          sparse_matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})>>(
          3., sparse::structural_zero, 4.)};

  assert((product.at<2, 1>() == 24.));
  assert((sum.at<0, 0>() == 3.));
  assert((difference.at<1, 0>() == -3.));
  assert((negated.at<0, 1>() == -2.));
  assert((scaled.at<0, 1>() == 4.));
  assert((divided.at<0, 1>() == 1.));
  assert((mirrored.at<0, 1>() == 3.));
  assert((widened == copy));
  assert(dot(column, column) == 25.);
  assert(magnitude(column) == 5.);
  assert(matrix_frob_norm(column) == 5.);
  assert(matrix_one_norm(band) == 12.);
  assert(matrix_inf_norm(band) == 13.);
  assert((band[1, 2] == 5.));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
