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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_EQUAL_TO_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_EQUAL_TO_TPP

namespace fcarouge {

//! @brief Element-wise equality of two typed matrices of the same shape.
//!
//! @details Compares each strongly typed element in turn, rather than the
//! backends' underlying storage: the storage holds bare representation
//! values, oblivious to the row and column indexes, so comparing it directly
//! would silently accept two matrices of incompatible element types (for
//! example, a matrix of lengths and a matrix of areas) whenever their raw
//! representations happen to match. Comparing elements routes back through
//! the strongly typed access, restoring the type check, and also sidesteps
//! backends, such as the `mdspan`-based one, whose underlying storage type
//! does not itself support equality.
//!
//! @note Deliberately excludes singleton, rank zero, matrices: those are
//! served by the dedicated overloads below.
[[nodiscard]] constexpr bool operator==(const rank_typed_matrix<2> auto &lhs,
                                        const rank_typed_matrix<2> auto &rhs) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;

  static_assert(same_shape<lhs_matrix, rhs_matrix>,
                "Matrix equality requires matrices of the same shapes, sizes.");

  bool result{true};

  tla::for_constexpr<lhs_matrix::rows>([&](auto i) {
    tla::for_constexpr<lhs_matrix::columns>([&](auto j) {
      using lhs_element = typename lhs_matrix::template element<i, j>;
      using rhs_element = typename rhs_matrix::template element<i, j>;

      static_assert(
          requires {
            std::declval<lhs_element>() == std::declval<rhs_element>();
          }, "Matrix equality requires comparable element types.");

      result &= (lhs.template at<i, j>() == rhs.template at<i, j>());
    });
  });

  return result;
}

//! @brief Element-wise equality of two typed matrices of the same shape.
//!
//! @see operator==(const rank_typed_matrix<2> auto &, const
//! rank_typed_matrix<2> auto &)
[[nodiscard]] constexpr bool operator==(const rank_typed_matrix<1> auto &lhs,
                                        const rank_typed_matrix<1> auto &rhs) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;

  static_assert(same_shape<lhs_matrix, rhs_matrix>,
                "Matrix equality requires matrices of the same shapes, sizes.");

  bool result{true};

  tla::for_constexpr<lhs_matrix::rows * lhs_matrix::columns>([&](auto i) {
    using lhs_element = typename lhs_matrix::template element<i>;
    using rhs_element = typename rhs_matrix::template element<i>;

    static_assert(
        requires {
          std::declval<lhs_element>() == std::declval<rhs_element>();
        }, "Matrix equality requires comparable element types.");

    result &= (lhs.template at<i>() == rhs.template at<i>());
  });

  return result;
}

//! @brief Equality of two singleton typed matrices.
[[nodiscard]] constexpr bool operator==(const rank_typed_matrix<0> auto &lhs,
                                        const rank_typed_matrix<0> auto &rhs) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;

  using lhs_element = typename lhs_matrix::template element<>;
  using rhs_element = typename rhs_matrix::template element<>;

  static_assert(
      requires { std::declval<lhs_element>() == std::declval<rhs_element>(); },
      "Matrix equality requires comparable element types.");

  return lhs.at() == rhs.at();
}

//! @brief Equality of a singleton typed matrix and another type.
//!
//! @details Symmetric equality, `lhs == rhs` and `rhs == lhs`, is
//! synthesized by the compiler from this single overload.
[[nodiscard]] constexpr bool operator==(const rank_typed_matrix<0> auto &lhs,
                                        const auto &rhs) {
  return lhs.at() == rhs;
}

} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_EQUAL_TO_TPP
