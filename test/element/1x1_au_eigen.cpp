/* Typed Linear Algebra
Version 0.3.0
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

#include <au/units/amperes.hh>
#include <au/units/meters.hh>
#include <au/units/newtons.hh>
#include <au/units/seconds.hh>

#include <cassert>
#include <concepts>
#include <tuple>

namespace fcarouge::test {
namespace {
//! @test Verifies the element type of the singleton, rank zero, matrix.
//!
//! @details The count of indexes must match the rank of the matrix: no
//! index for a singleton matrix.
[[maybe_unused]] const auto test{[] -> int {
  using matrix = matrix<double, std::tuple<decltype(au::meters(1.))>,
                        std::tuple<decltype(au::meters(1.))>>;
  using quantity = decltype(au::squared(au::meters)(1.));

  static_assert(std::same_as<matrix::element<>, quantity>);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
