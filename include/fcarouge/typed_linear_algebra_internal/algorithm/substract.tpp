/* Typed Linear Algebra
Version 0.2.0
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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_SUBSTRACT_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_SUBSTRACT_TPP

namespace fcarouge {

//! @todo Requires, assert that the element types are compatible.
[[nodiscard]] constexpr auto operator-(const same_as_typed_matrix auto &lhs,
                                       const same_as_typed_matrix auto &rhs) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;

  static_assert(
      same_shape<lhs_matrix, rhs_matrix>,
      "Matrix subtraction requires matrices of the same shapes, sizes.");

  // Each typed element of the lhs matrix must be substractable to the
  // corresponding typed element of the rhs matrix.
  tla::for_constexpr<0, lhs_matrix::rows, 1>([&](auto i) {
    tla::for_constexpr<0, lhs_matrix::columns, 1>([&](auto j) {
      using lhs_element = typename lhs_matrix::template element<i, j>;
      using rhs_element = typename rhs_matrix::template element<i, j>;

      static_assert(
          requires {
            std::declval<lhs_element>() - std::declval<rhs_element>();
          }, "Matrix subtraction requires compatible element types.");
    });
  });

  using row_indexes = typename lhs_matrix::row_indexes;
  using column_indexes = typename lhs_matrix::column_indexes;

  if constexpr (lhs_matrix::rank > 0) {
    return make_typed_matrix<row_indexes, column_indexes>(lhs.data() -
                                                          rhs.data());
  } else {
    using lhs_element = typename lhs_matrix::template element<0, 0>;
    using rhs_element = typename rhs_matrix::template element<0, 0>;

    return lhs_element{lhs} - rhs_element{rhs};
  }
}

[[nodiscard]] constexpr auto operator-(const other auto &lhs,
                                       const singleton_typed_matrix auto &rhs) {
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<0, 0>;

  return lhs - element{rhs};
}

[[nodiscard]] constexpr auto operator-(const singleton_typed_matrix auto &lhs,
                                       const other auto &rhs) {
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<0, 0>;

  return element{lhs} - rhs;
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_SUBSTRACT_TPP
