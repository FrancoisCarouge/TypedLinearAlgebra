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
#include <utility>

namespace fcarouge::test {
namespace {
//! @test Verifies the construction from an underlying matrix of every value
//! category.
[[maybe_unused]] const auto test{[] -> int {
  using typed = matrix<double, 2, 3>;
  using storage = typed::matrix;

  const typed expected{{1., 2., 3.}, {4., -5., 6.}};

  storage lvalue;
  const storage const_lvalue{expected.data()};
  storage xvalue{expected.data()};

  lvalue = expected.data();

  const typed from_lvalue{lvalue};
  const typed from_const_lvalue{const_lvalue};
  const typed from_xvalue{std::move(xvalue)};
  const typed from_prvalue{storage{expected.data()}};

  assert(from_lvalue == expected);
  assert(from_const_lvalue == expected);
  assert(from_xvalue == expected);
  assert(from_prvalue == expected);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
