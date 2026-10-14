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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_RANK_1_UPDATE_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_RANK_1_UPDATE_TPP

//! @todo Remove the feature check when supporting native C++26.
#ifdef __cpp_lib_linalg

#include "fcarouge/typed_linear_algebra/internal/as_vector_span.hpp"

#include <linalg>

namespace fcarouge {

//! @brief Overwriting rank-1 update, the outer product `result = lhs rhsᵀ` of
//! two vectors.
//!
//! @details The element of the result at row `i` and column `j` is the
//! product of the `i`-th element of `lhs` and the `j`-th element of `rhs`, and
//! must be assignable to the corresponding typed element of `result`. A matrix
//! whose row indexes are those of `lhs` and whose column indexes are those of
//! `rhs` is the natural result type. Row and column vectors are both accepted
//! for `lhs` and `rhs`, as only their elements' order matters.
//!
//! @see std::linalg::matrix_rank_1_update
constexpr void matrix_rank_1_update(const rank_typed_matrix<1> auto &lhs,
                                    const rank_typed_matrix<1> auto &rhs,
                                    rank_typed_matrix<2> auto &result) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;
  using result_matrix = std::remove_cvref_t<decltype(result)>;

  static_assert(
      result_matrix::rows == lhs_matrix::rows * lhs_matrix::columns,
      "Matrix rank-1 update requires the lhs vector size to match the "
      "matrix row count.");
  static_assert(
      result_matrix::columns == rhs_matrix::rows * rhs_matrix::columns,
      "Matrix rank-1 update requires the rhs vector size to match the "
      "matrix column count.");

  // Each product of the typed elements of the lhs and rhs vectors must be
  // assignable to the corresponding typed element of the result matrix.
  tla::for_constexpr<result_matrix::rows>([&](auto i) {
    tla::for_constexpr<result_matrix::columns>([&](auto j) {
      using lhs_element = typename lhs_matrix::template element<i>;
      using rhs_element = typename rhs_matrix::template element<j>;
      using result_element = typename result_matrix::template element<i, j>;

      static_assert(
          requires {
            std::declval<result_element &>() =
                std::declval<lhs_element>() * std::declval<rhs_element>();
          }, "Matrix rank-1 update requires compatible element types.");
    });
  });

  using std::linalg::matrix_rank_1_update;
  matrix_rank_1_update(tla::as_vector_span(lhs), tla::as_vector_span(rhs),
                       result.data());
}

//! @brief Updating rank-1 update, the sum `result = addend + lhs rhsᵀ` of a
//! matrix and the outer product of two vectors.
//!
//! @details The element of the result at row `i` and column `j` is the sum
//! of the corresponding typed element of `addend` and the product of the `i`-th
//! element of `lhs` and the `j`-th element of `rhs`, and must be assignable to
//! the corresponding typed element of `result`. The `addend` and `result`
//! matrices may be the same object, for an in-place update
//! `result += lhs rhsᵀ`.
//!
//! @see std::linalg::matrix_rank_1_update
constexpr void matrix_rank_1_update(const rank_typed_matrix<1> auto &lhs,
                                    const rank_typed_matrix<1> auto &rhs,
                                    const rank_typed_matrix<2> auto &addend,
                                    rank_typed_matrix<2> auto &result) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;
  using addend_matrix = std::remove_cvref_t<decltype(addend)>;
  using result_matrix = std::remove_cvref_t<decltype(result)>;

  static_assert(
      result_matrix::rows == lhs_matrix::rows * lhs_matrix::columns,
      "Matrix rank-1 update requires the lhs vector size to match the "
      "matrix row count.");
  static_assert(
      result_matrix::columns == rhs_matrix::rows * rhs_matrix::columns,
      "Matrix rank-1 update requires the rhs vector size to match the "
      "matrix column count.");
  static_assert(
      same_shape<addend_matrix, result_matrix>,
      "Matrix rank-1 update requires the addend and result matrices of the "
      "same shapes, sizes.");

  // Each typed element of the addend matrix must be addable to the product of
  // the corresponding typed elements of the lhs and rhs vectors, and the sum
  // must be assignable to the corresponding typed element of the result matrix.
  tla::for_constexpr<result_matrix::rows>([&](auto i) {
    tla::for_constexpr<result_matrix::columns>([&](auto j) {
      using lhs_element = typename lhs_matrix::template element<i>;
      using rhs_element = typename rhs_matrix::template element<j>;
      using addend_element = typename addend_matrix::template element<i, j>;
      using result_element = typename result_matrix::template element<i, j>;

      static_assert(
          requires {
            std::declval<result_element &>() =
                std::declval<addend_element>() +
                std::declval<lhs_element>() * std::declval<rhs_element>();
          }, "Matrix rank-1 update requires compatible element types.");
    });
  });

  using std::linalg::matrix_rank_1_update;
  matrix_rank_1_update(tla::as_vector_span(lhs), tla::as_vector_span(rhs),
                       addend.data(), result.data());
}
} // namespace fcarouge

#endif
#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_RANK_1_UPDATE_TPP
