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
#include <cstddef>
#include <utility>

namespace fcarouge::test {
namespace {
//! @test Verifies the operations on a large typed matrix, on run-time values,
//! in a real-time context: no memory allocation, no lock, no blocking call.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  constexpr std::size_t size{32};

  using band_matrix = sparse_matrix<double, sparse::tridiagonal<size>>;

  band_matrix band;

  // A one-dimensional Laplacian.
  [&band]<std::size_t... Indexes>(std::index_sequence<Indexes...>) -> void {
    (band.at<Indexes, Indexes>(2.), ...);
    (band.at<Indexes, Indexes + 1>(-1.), ...);
    (band.at<Indexes + 1, Indexes>(-1.), ...);
  }(std::make_index_sequence<size - 1>{});
  band.at<size - 1, size - 1>(2.);

  const auto square{band * band};

  // Far from the boundaries: the stencil squared.
  assert((square.at<10, 10>() == 6.));
  assert((square.at<10, 11>() == -4.));
  assert((square.at<10, 12>() == 1.));
  assert((square.at<10, 13>() == 0.));
  assert((square.at<0, 0>() == 5.));
  assert((square.at<31, 31>() == 5.));

  // A diagonal scaling, a sum, a transposition, and the norms.
  const auto sum{band + transposed(band)};

  assert((sum.at<10, 10>() == 4.));
  assert((sum.at<10, 11>() == -2.));
  assert((matrix_inf_norm(band) == 4.));
  assert((matrix_one_norm(band) == 4.));

  // A matrix-vector product: the Laplacian of a constant vanishes inside.
  sparse_matrix<double, sparse::dense<size, 1>> constant;

  [&constant]<std::size_t... Indexes>(std::index_sequence<Indexes...>) -> void {
    (constant.at<Indexes>(1.), ...);
  }(std::make_index_sequence<size>{});

  const auto image{band * constant};

  assert(image.at<0>() == 1.);
  assert(image.at<10>() == 0.);
  assert(image.at<31>() == 1.);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
