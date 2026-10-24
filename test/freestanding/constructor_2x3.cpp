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

namespace fcarouge::test {
namespace {
//! @test Verifies the construction and the element access of a matrix on a
//! freestanding implementation.
static_assert([] -> bool {
  const matrix<double, 2, 3> zero{};
  matrix<double, 2, 3> value{{1., 2., 3.}, {4., 5., 6.}};

  value.at<1, 2>(-6.);

  return zero.at<0, 0>() == 0. && zero.at<1, 2>() == 0. &&
         value.at<0, 0>() == 1. && value.at<0, 1>() == 2. &&
         value.at<0, 2>() == 3. && value.at<1, 0>() == 4. &&
         value.at<1, 1>() == 5. && value.at<1, 2>() == -6.;
}());

//! @test Verifies the construction and the element access of a vector on a
//! freestanding implementation.
static_assert([] -> bool {
  const column_vector<double, 3> value{1., 2., 3.};

  return value.at<0>() == 1. && value.at<1>() == 2. && value.at<2>() == 3.;
}());
} // namespace
} // namespace fcarouge::test
