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
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::one;

namespace {
//! @test Verifies the multiplication operators of a singleton matrix of a
//! dimensionless element, whose quantity is also constructible from a number.
[[maybe_unused]] const auto test{[] -> int {
  using ratio = quantity<mp_units::dimensionless[one]>;

  const ratio lhs_value{2. * one};
  const ratio rhs_value{3. * one};
  const ratio scalar{4. * one};

  const row_vector<representation, ratio> lhs{lhs_value};
  const row_vector<representation, ratio> rhs{rhs_value};

  const ratio both{lhs * rhs};
  const ratio scalar_lhs{scalar * rhs};
  const ratio scalar_rhs{lhs * scalar};

  assert(6. * one == both);
  assert(12. * one == scalar_lhs);
  assert(8. * one == scalar_rhs);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
