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
#include <chrono>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the dot product of a row vector of std::chrono durations
//! and a row vector of the dimensionless representation type: durations have
//! no duration-by-duration product, so each term instead pairs a duration
//! with a plain scalar, mirroring `matrix_vector_product`'s convention.
[[maybe_unused]] const auto test{[] -> int {
  using seconds = std::chrono::duration<representation>;

  const row_vector<representation, seconds, seconds> a{seconds{2.},
                                                       seconds{3.}};
  const row_vector<representation, representation, representation> b{4., 5.};

  assert(dot(a, b) == seconds{23.});

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
