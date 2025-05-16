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
  using type = decltype(value);
  using matrix = std::remove_cvref_t<type>;

  return typed_matrix<matrix, RowIndexes, ColumnIndexes>{
      std::forward<type>(value)};
}

//! @todo Requires, assert that the element types are compatible.
[[nodiscard]] constexpr bool operator==(const is_typed_matrix auto &lhs,
                                        const is_typed_matrix auto &rhs) {
  return lhs.data() == rhs.data();
}

template <typename Matrix1, typename RowIndexes1, typename ColumnIndexes1,
          typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
[[nodiscard]] constexpr auto
operator*(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &rhs) {
  //! @todo Simplify the size check with tla::size.
  static_assert(typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1>::columns ==
                    typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2>::rows,
                "Matrix multiplication requires compatible sizes.");

  // A valid combination of row, column indexes providing the expected
  // element types. Any Nth combination, here 0th, of indexes is expected to
  // be valid.
  using row_indexes =
      tla::product<RowIndexes1, std::tuple_element_t<0, ColumnIndexes1>>;
  using column_indexes =
      tla::product<ColumnIndexes2, std::tuple_element_t<0, RowIndexes2>>;

  // Verify every terms of the inner product for every elements is of the same
  // type. For each of the column indexes of the left-hand-side matrix factor.
  tla::for_constexpr<0, std::tuple_size_v<ColumnIndexes1>, 1>([&](auto i) {
    // Find the type of the first factor of the Ith term of the element's inner
    // product. That is, the Ith row index tuple.
    using FirstFactorIthSum =
        tla::product<RowIndexes1, std::tuple_element_t<i, ColumnIndexes1>>;
    // For each of the row indexes of the right-hand-side matrix factor.
    tla::for_constexpr<0, std::tuple_size_v<RowIndexes2>, 1>([&](auto j) {
      // Find the type of the second factor of the Jth term of the element's
      // inner product. That is, the Jth column index tuple.
      using SecondFactorJthSum =
          tla::product<ColumnIndexes2, std::tuple_element_t<j, RowIndexes2>>;
      // For each term of the inner product of the element.
      tla::for_constexpr<0, std::tuple_size_v<RowIndexes1>, 1>([&](auto k) {
        tla::for_constexpr<0, std::tuple_size_v<ColumnIndexes2>, 1>(
            [&](auto l) {
              static_assert(
                  std::is_same_v<
                      tla::product<std::tuple_element_t<k, FirstFactorIthSum>,
                                   std::tuple_element_t<l, SecondFactorJthSum>>,
                      tla::product<std::tuple_element_t<k, row_indexes>,
                                   std::tuple_element_t<l, column_indexes>>>,
                  "Matrix multiplication requires compatible types.");
            });
      });
    });
  });

  return make_typed_matrix<row_indexes, column_indexes>(lhs.data() *
                                                        rhs.data());
}

[[nodiscard]] constexpr auto operator*(const is_typed_matrix auto &lhs,
                                       const auto &rhs) {
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(rhs)>;
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using row_indexes = typename matrix::row_indexes;
  using column_indexes = typename matrix::column_indexes;
  using underlying = typename matrix::underlying;

  return make_typed_matrix<tla::product<row_indexes, type>, column_indexes>(
      lhs.data() * cast<underlying, type>(rhs));
}

[[nodiscard]] constexpr auto operator*(const auto &lhs,
                                       const is_typed_matrix auto &rhs) {
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(lhs)>;
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using row_indexes = typename matrix::row_indexes;
  using column_indexes = typename matrix::column_indexes;
  using underlying = typename matrix::underlying;

  return make_typed_matrix<tla::product<row_indexes, type>, column_indexes>(
      cast<underlying, type>(lhs) * rhs.data());
}

[[nodiscard]] constexpr auto operator*(const auto &lhs,
                                       const is_typed_matrix auto &rhs)
  requires is_singleton_typed_matrix<decltype(rhs)>
{
  //! @todo Should there be constraints on the type?
  using matrix = decltype(rhs);
  using element = typename matrix::template element<0, 0>;

  return lhs * element{rhs};
}

[[nodiscard]] constexpr auto operator*(const is_typed_matrix auto &lhs,
                                       const auto &rhs)
  requires is_singleton_typed_matrix<decltype(lhs)>
{
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using element = typename matrix::template element<0, 0>;

  return element{lhs} * rhs;
}

//! @todo Requires, assert that the element types are compatible.
[[nodiscard]] constexpr auto operator+(const is_typed_matrix auto &lhs,
                                       const is_typed_matrix auto &rhs) {
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using row_indexes = typename matrix::row_indexes;
  using column_indexes = typename matrix::column_indexes;

  return make_typed_matrix<row_indexes, column_indexes>(lhs.data() +
                                                        rhs.data());
}

[[nodiscard]] constexpr auto operator+(const auto &lhs,
                                       const is_typed_matrix auto &rhs)
  requires is_singleton_typed_matrix<decltype(rhs)>
{
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<0, 0>;

  return lhs + element{rhs};
}

[[nodiscard]] constexpr auto operator+(const is_typed_matrix auto &lhs,
                                       const auto &rhs)
  requires is_singleton_typed_matrix<decltype(lhs)>
{
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<0, 0>;

  return element{lhs} + rhs;
}

//! @todo Requires, assert that the element types are compatible.
[[nodiscard]] constexpr auto operator-(const is_typed_matrix auto &lhs,
                                       const is_typed_matrix auto &rhs) {
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using row_indexes = typename matrix::row_indexes;
  using column_indexes = typename matrix::column_indexes;

  return make_typed_matrix<row_indexes, column_indexes>(lhs.data() -
                                                        rhs.data());
}

[[nodiscard]] constexpr auto operator-(const auto &lhs,
                                       const is_typed_matrix auto &rhs)
  requires is_singleton_typed_matrix<decltype(rhs)>
{
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<0, 0>;

  return lhs - element{rhs};
}

[[nodiscard]] constexpr auto operator-(const is_typed_matrix auto &lhs,
                                       const auto &rhs)
  requires is_singleton_typed_matrix<decltype(lhs)>
{
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<0, 0>;

  return element{lhs} - rhs;
}

template <typename Matrix1, typename RowIndexes1, typename ColumnIndexes1,
          typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
[[nodiscard]] constexpr auto
operator/(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &rhs)
  requires is_column_typed_matrix<
      typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2>>
{
  //! @todo Simplify the size check with tla::size.
  static_assert(typed_matrix<Matrix1, RowIndexes1, ColumnIndexes1>::rows ==
                    typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2>::columns,
                "Matrix division requires compatible sizes.");

  //! @todo This is likely incorrect, incomplete?
  using row_indexes =
      tla::quotient<std::tuple_element_t<0, RowIndexes1>, ColumnIndexes2>;
  using column_indexes =
      tla::quotient<std::tuple_element_t<0, ColumnIndexes1>, RowIndexes2>;

  //! @todo Add type verification, perhaps with a generalization of the
  //! multiplication verification?

  return make_typed_matrix<row_indexes, column_indexes>(lhs.data() /
                                                        rhs.data());
}

[[nodiscard]] constexpr auto operator/(const is_typed_matrix auto &lhs,
                                       const auto &rhs) {
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(rhs)>;
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using underlying = typename matrix::underlying;
  using row_indexes = tla::quotient<typename matrix::row_indexes, type>;
  using column_indexes = typename matrix::column_indexes;

  return make_typed_matrix<row_indexes, column_indexes>(
      lhs.data() / cast<underlying, type>(rhs));
}

[[nodiscard]] constexpr auto operator/(const is_typed_matrix auto &lhs,
                                       const auto &rhs)
  requires is_singleton_typed_matrix<decltype(lhs)>
{
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using element = typename matrix::template element<0, 0>;

  return element{lhs} / rhs;
}

[[nodiscard]] constexpr auto operator/(const auto &lhs,
                                       const is_typed_matrix auto &rhs)
  requires is_column_typed_matrix<decltype(rhs)>
{
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(lhs)>;
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using underlying = typename matrix::underlying;
  using row_indexes = tla::quotient<type, typename matrix::column_indexes>;
  using column_indexes =
      tla::quotient<std::identity, typename matrix::row_indexes>;

  //! @todo Add type verification, perhaps with a generalization of the
  //! multiplication verification?

  return make_typed_matrix<row_indexes, column_indexes>(
      cast<underlying, type>(lhs) / rhs.data());
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_OPERATION_TPP
