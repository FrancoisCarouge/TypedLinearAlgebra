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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_PRODUCT_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_PRODUCT_TPP

//! @todo Remove the feature check when supporting native C++26.
#ifdef __cpp_lib_linalg

#include <linalg>

namespace fcarouge {

//! @brief Computes the product of two matrices.
//!
//! @details Every term summed into an element of the result must be of the
//! scale of that element: the stored numbers are computed on directly. Convert
//! the operands first otherwise.
//!
//! @see std::linalg::matrix_product
constexpr void matrix_product(const same_as_typed_matrix auto &lhs,
                              const same_as_typed_matrix auto &rhs,
                              same_as_typed_matrix auto &result) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;
  using result_matrix = std::remove_cvref_t<decltype(result)>;

  tla::for_constexpr<lhs_matrix::rows>([](auto i) {
    tla::for_constexpr<rhs_matrix::columns>([i](auto j) {
      tla::for_constexpr<lhs_matrix::columns>([i, j](auto k) {
        static_cast<void>(i); // Compiler compatibility.
        static_cast<void>(j); // Compiler compatibility.
        static_assert(
            tla::is_same_scale<tla::product<tla::element_at<lhs_matrix, i, k>,
                                            tla::element_at<rhs_matrix, k, j>>,
                               tla::element_at<result_matrix, i, j>,
                               typename lhs_matrix::underlying>(),
            "Matrix product requires a result of the scale of its terms, "
            "convert first.");
      });
    });
  });

  using std::linalg::matrix_product;
  matrix_product(lhs.data(), rhs.data(), result.data());
}
} // namespace fcarouge

#endif
#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_PRODUCT_TPP
