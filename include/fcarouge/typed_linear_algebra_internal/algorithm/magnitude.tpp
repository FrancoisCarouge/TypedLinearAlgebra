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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_MAGNITUDE_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_MAGNITUDE_TPP

#include <cmath>

namespace fcarouge {
[[nodiscard]] constexpr auto magnitude(const uniform_typed_matrix auto &value) {
  // There exists definitions of magnitude to support for other shapes/ranks.
  static_assert(
      rank_typed_matrix<decltype(value), 1>,
      "The magnitude operation only supports vector types at this time.");

  using matrix = std::remove_cvref_t<decltype(value)>;
  using element = typename matrix::template element<0>;
  using underlying = typename matrix::underlying;

  underlying sums{};

  // There exists a variety of implementation tradeoffs to explore. Delegate to
  // underlying linear algebra library? Implement atop strong types?
  tla::for_constexpr<matrix::rows>([&](auto i) {
    sums += cast<underlying, element>(value.template at<i>()) *
            cast<underlying, element>(value.template at<i>());
  });

  using std::sqrt;

  return cast<element, underlying>(sqrt(sums));
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_MAGNITUDE_TPP
