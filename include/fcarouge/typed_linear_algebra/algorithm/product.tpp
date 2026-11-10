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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_PRODUCT_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_PRODUCT_TPP

namespace fcarouge {
//! @brief Concept of typed matrix shapes admitting a matrix product.
//!
//! @details Excludes singleton by singleton, an ordinary scalar product
//! served by a dedicated, simpler overload.
template <typename Lhs, typename Rhs>
concept multipliable_shape =
    same_as_typed_matrix<Lhs> and same_as_typed_matrix<Rhs> and
    (std::remove_cvref_t<Lhs>::columns == std::remove_cvref_t<Rhs>::rows) and
    (std::remove_cvref_t<Lhs>::columns > 1 or
     std::remove_cvref_t<Lhs>::rows > 1 or
     std::remove_cvref_t<Rhs>::columns > 1);

//! @brief Whether every per-term product of the `Lhs * Rhs` matrix product
//! converts to the type of its row-column's first term, as required to sum
//! them.
//!
//! @details The terms of the first row-column are verified, at a cost linear in
//! the inner dimension rather than cubic in the size of the matrices: an
//! element type is the product of its row and column indexes, so each term of
//! another row-column is the corresponding term of the first row-column times
//! the same row and column index factors. This holds for element types whose
//! product is associative and commutative, as the result indexes, formed from
//! the first column and row of the operands, already require. Matrices of
//! `std::identity` indexes only have no terms to verify.
//!
//! @pre `Lhs` and `Rhs` are `multipliable_shape`.
template <typename Lhs, typename Rhs>
constexpr auto are_terms_multipliable() -> bool {
  using lhs_matrix = std::remove_cvref_t<Lhs>;
  using rhs_matrix = std::remove_cvref_t<Rhs>;

  if constexpr (tla::identity_typed_matrix<lhs_matrix> and
                tla::identity_typed_matrix<rhs_matrix>) {
    return true;
  } else {
    using first = tla::product<tla::element_at<lhs_matrix, 0, 0>,
                               tla::element_at<rhs_matrix, 0, 0>>;

    bool convertible{true};

    tla::for_constexpr<lhs_matrix::columns>([&convertible](auto k) {
      convertible &=
          std::is_convertible_v<tla::product<tla::element_at<lhs_matrix, 0, k>,
                                             tla::element_at<rhs_matrix, k, 0>>,
                                first>;
    });

    return convertible;
  }
}

//! @brief Concept of typed matrices whose per-term products can be summed
//! into a matrix product result.
template <typename Lhs, typename Rhs>
concept multipliable_elements = are_terms_multipliable<Lhs, Rhs>();

//! @brief Whether every per-term product of the `Lhs * Rhs` matrix product is
//! of the scale of its row-column's first term, as required to sum their stored
//! numbers.
//!
//! @details The terms of the first row-column are verified, for the reasons
//! of `are_terms_multipliable`: the scale of a product is the product of the
//! scales of its factors.
//!
//! @pre `Lhs` and `Rhs` are `multipliable_shape`.
template <typename Lhs, typename Rhs>
constexpr auto are_terms_same_scale() -> bool {
  using lhs_matrix = std::remove_cvref_t<Lhs>;
  using rhs_matrix = std::remove_cvref_t<Rhs>;

  if constexpr (tla::identity_typed_matrix<lhs_matrix> and
                tla::identity_typed_matrix<rhs_matrix>) {
    return true;
  } else {
    using underlying = typename lhs_matrix::underlying;
    using first = tla::product<tla::element_at<lhs_matrix, 0, 0>,
                               tla::element_at<rhs_matrix, 0, 0>>;

    bool same_scale{true};

    tla::for_constexpr<lhs_matrix::columns>([&same_scale](auto k) {
      same_scale &=
          tla::is_same_scale<tla::product<tla::element_at<lhs_matrix, 0, k>,
                                          tla::element_at<rhs_matrix, k, 0>>,
                             first, underlying>();
    });

    return same_scale;
  }
}

//! @brief Concept of typed matrices whose per-term products can be summed from
//! their stored numbers: convert the operands first otherwise.
template <typename Lhs, typename Rhs>
concept same_scale_terms = are_terms_same_scale<Lhs, Rhs>();

//! @brief Concept of typed matrices that can be multiplied together.
template <typename Lhs, typename Rhs>
concept multipliable =
    multipliable_shape<Lhs, Rhs> and multipliable_elements<Lhs, Rhs> and
    same_scale_terms<Lhs, Rhs>;

//! @brief Concept of a type scalable by another.
//!
//! @details Named so the deduced-return-type overloads below can constrain
//! themselves on it, rather than leaving the multiplication unconstrained: an
//! unconstrained, deduced return type forces the compiler to instantiate the
//! function body, and thus fail hard, merely to probe whether the
//! multiplication is well-formed, for example when third party code, such as
//! `mp-units`, uses a `requires` expression to check whether some unrelated
//! type supports multiplication by a rank zero typed matrix's element type.
template <typename Scaled, typename Scale>
concept scalable_by = requires {
  std::declval<const Scaled &>() * std::declval<const Scale &>();
};

[[nodiscard]] constexpr auto operator*(const same_as_typed_matrix auto &lhs,
                                       const same_as_typed_matrix auto &rhs)
  requires multipliable<decltype(lhs), decltype(rhs)>
{
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;
  using lhs_column_indexes = typename lhs_matrix::column_indexes;
  using rhs_row_indexes = typename rhs_matrix::row_indexes;
  using row_indexes = tla::product<typename lhs_matrix::row_indexes,
                                   std::tuple_element_t<0, lhs_column_indexes>>;
  using column_indexes = tla::product<typename rhs_matrix::column_indexes,
                                      std::tuple_element_t<0, rhs_row_indexes>>;

  return make_typed_matrix<row_indexes, column_indexes>(lhs.data() *
                                                        rhs.data());
}

[[nodiscard]] constexpr auto operator*(const same_as_typed_matrix auto &lhs,
                                       const other auto &rhs) {
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(rhs)>;
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using row_indexes = typename matrix::row_indexes;
  using column_indexes = typename matrix::column_indexes;
  using underlying = typename matrix::underlying;

  return make_typed_matrix<tla::product<row_indexes, type>, column_indexes>(
      lhs.data() * cast<underlying, type>(rhs));
}

[[nodiscard]] constexpr auto operator*(const other auto &lhs,
                                       const same_as_typed_matrix auto &rhs) {
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(lhs)>;
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using row_indexes = typename matrix::row_indexes;
  using column_indexes = typename matrix::column_indexes;
  using underlying = typename matrix::underlying;

  return make_typed_matrix<tla::product<row_indexes, type>, column_indexes>(
      cast<underlying, type>(lhs) * rhs.data());
}

[[nodiscard]] constexpr auto operator*(const other auto &lhs,
                                       const rank_typed_matrix<0> auto &rhs)
  requires scalable_by<decltype(lhs), typename std::remove_cvref_t<
                                          decltype(rhs)>::template element<>>
{
  return lhs * rhs.template at<>();
}

//! @see operator*(const other auto &, const rank_typed_matrix<0> auto &)
[[nodiscard]] constexpr auto operator*(const rank_typed_matrix<0> auto &lhs,
                                       const other auto &rhs)
  requires scalable_by<
      typename std::remove_cvref_t<decltype(lhs)>::template element<>,
      decltype(rhs)>
{
  return lhs.template at<>() * rhs;
}

[[nodiscard]] constexpr auto operator*(const rank_typed_matrix<0> auto &lhs,
                                       const rank_typed_matrix<0> auto &rhs) {
  return lhs.template at<>() * rhs.template at<>();
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_PRODUCT_TPP
