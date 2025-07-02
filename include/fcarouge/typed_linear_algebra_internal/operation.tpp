/* Typed Linear Algebra
Version 0.1.0
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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_OPERATION_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_OPERATION_TPP

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

//! @brief Factory function in support of syntaxic reduction.
//!
//! @todo Should this be an internal only tool? Should, can there be a safe
//! public version?
//!
//! @warning Useful for operations implementation where underlying data
//! constrution is needed. Not recommended for convenience construction due to
//! absence of type validation.
template <typename RowIndexes, typename ColumnIndexes>
auto make_typed_matrix(auto &&value) {
  return typed_matrix<std::remove_cvref_t<decltype(value)>, RowIndexes,
                      ColumnIndexes>{std::forward<decltype(value)>(value)};
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr bool
operator==(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
           const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return lhs.data() == rhs.data();
}

template <typename Matrix1, typename RowIndexes1, typename ColumnIndexes1,
          typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &rhs) {

  //! @todo Simplify the size check with tla::size.
  static_assert(typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1>::columns ==
                    typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2>::rows,
                "Matrix multiplication requires compatible sizes.");

  // A valid combination of row, column indexes providing the expected
  // element types. Any Nth combination, here 0th, of indexes is expected to
  // be valid.
  using RowIndexes =
      tla::product<RowIndexes1, std::tuple_element_t<0, ColumnIndexes1>>;
  using ColumnIndexes =
      tla::product<ColumnIndexes2, std::tuple_element_t<0, RowIndexes2>>;
  /*
  // Verify every terms of the inner product for every elements is of the same
  // type. For each of the column indexes of the left-hand-side matrix factor.
  tla::for_constexpr<0, tla::size<ColumnIndexes1>, 1>([&](auto i) {
    // Find the type of the first factor of the Ith term of the element's inner
    // product. That is, the Ith row index tuple.
    using FirstFactorIthSum =
        tla::product<RowIndexes1, std::tuple_element_t<i, ColumnIndexes1>>;
    // For each of the row indexes of the right-hand-side matrix factor.
    tla::for_constexpr<0, tla::size<RowIndexes2>, 1>([&](auto j) {
      // Find the type of the second factor of the Jth term of the element's
      // inner product. That is, the Jth column index tuple.
      using SecondFactorJthSum =
          tla::product<ColumnIndexes2, std::tuple_element_t<j, RowIndexes2>>;
      // For each term of the inner product of the element.
      tla::for_constexpr<0, tla::size<RowIndexes1>, 1>([&](auto k) {
        tla::for_constexpr<0, tla::size<ColumnIndexes2>, 1>([&](auto l) {
          static_assert(
              std::is_same_v<
                  tla::product<std::tuple_element_t<k, FirstFactorIthSum>,
                               std::tuple_element_t<l, SecondFactorJthSum>>,
                  tla::product<std::tuple_element_t<k, RowIndexes>,
                               std::tuple_element_t<l, ColumnIndexes>>>,
              "Matrix multiplication requires compatible types.");
        });
      });
    });
  });
  */
  return make_typed_matrix<RowIndexes, ColumnIndexes>(lhs.data() * rhs.data());
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes,
          typename Type>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          const Type &rhs) {
  using matrix = typed_matrix<Matrix, RowIndexes, ColumnIndexes>;
  using underlying = typename matrix::underlying;
  return make_typed_matrix<tla::product<RowIndexes, Type>, ColumnIndexes>(
      lhs.data() * cast<underlying, Type>(rhs));
}

template <typename Type, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator*(const Type &lhs,
          const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &rhs) {
  using matrix = typed_matrix<Matrix, RowIndexes, ColumnIndexes>;
  using underlying = typename matrix::underlying;
  return make_typed_matrix<tla::product<RowIndexes, Type>, ColumnIndexes>(
      cast<underlying, Type>(lhs) * rhs.data());
}

template <typename Type, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
[[nodiscard]] inline constexpr auto
operator*(const Type &lhs,
          const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &rhs) {
  return lhs *
         tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, 0, 0>{
             rhs.data()};
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes,
          typename Type>
  requires tla::singleton<typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          const Type &rhs) {
  return tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, 0, 0>{
             lhs.data()} *
         rhs;
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator+(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return make_typed_matrix<RowIndexes, ColumnIndexes>(lhs.data() + rhs.data());
}

template <typename Matrix1, typename RowIndexes1, typename ColumnIndexes1,
          typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
[[nodiscard]] inline constexpr auto
operator+(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &rhs) {
  //! @todo Verify the resulting types compatibility at compile-time.
  return make_typed_matrix<RowIndexes1, ColumnIndexes1>(lhs.data() +
                                                        rhs.data());
}

//! @todo Generalize out the scalar restriction.
template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
[[nodiscard]] inline constexpr auto
operator+(Scalar lhs,
          const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &rhs) {
  return lhs +
         tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, 0, 0>{
             rhs.data()};
}

//! @todo Generalize out the scalar restriction.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes,
          tla::arithmetic Scalar>
  requires tla::singleton<typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
[[nodiscard]] inline constexpr auto
operator+(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, 0, 0>{
             lhs.data()} +
         rhs;
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator-(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return make_typed_matrix<RowIndexes, ColumnIndexes>(lhs.data() - rhs.data());
}

template <typename Matrix1, typename RowIndexes1, typename ColumnIndexes1,
          typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
[[nodiscard]] inline constexpr auto
operator-(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &rhs) {
  //! @todo Verify the resulting types compatibility at compile-time.
  return make_typed_matrix<RowIndexes1, ColumnIndexes1>(lhs.data() -
                                                        rhs.data());
}

//! @todo Generalize out the scalar restriction.
template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
[[nodiscard]] inline constexpr auto
operator-(Scalar lhs,
          const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &rhs) {
  return lhs -
         tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, 0, 0>{
             rhs.data()};
}

//! @todo Generalize out the scalar restriction.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes,
          tla::arithmetic Scalar>
  requires tla::singleton<typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
[[nodiscard]] inline constexpr auto
operator-(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, 0, 0>{
             lhs.data()} -
         rhs;
}

template <typename Matrix1, typename Matrix2, typename RowIndexes1,
          typename ColumnIndexes1, typename RowIndexes2,
          typename ColumnIndexes2>
[[nodiscard]] inline constexpr auto
operator/(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &rhs) {

  //! @todo Simplify the size check with tla::size.
  static_assert(typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1>::columns ==
                    typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2>::columns,
                "Matrix division requires compatible sizes.");

  //! @todo This is likely incorrect, incomplete?
  using RowIndexes =
      tla::quotient<RowIndexes1, std::tuple_element_t<0, RowIndexes2>>;
  using ColumnIndexes =
      tla::quotient<std::tuple_element_t<0, RowIndexes1>, RowIndexes2>;

  //! @todo Add type verification, perhaps with a generalization of the
  //! multiplication verification?

  return make_typed_matrix<RowIndexes, ColumnIndexes>(lhs.data() / rhs.data());
}

//! @todo Generalize out the scalar restriction.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes,
          tla::arithmetic Scalar>
[[nodiscard]] inline constexpr auto
operator/(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return make_typed_matrix<RowIndexes, ColumnIndexes>(lhs.data() / rhs);
}

//! @todo Generalize out the scalar restriction.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes,
          tla::arithmetic Scalar>
  requires tla::singleton<typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
[[nodiscard]] inline constexpr auto
operator/(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, 0, 0>{
             lhs.data()} /
         rhs;
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_OPERATION_TPP
