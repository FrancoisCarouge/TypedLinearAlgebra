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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_PRODUCT_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_PRODUCT_TPP

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

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
//! @pre `Lhs` and `Rhs` are `multipliable_shape`.
template <typename Lhs, typename Rhs> constexpr bool are_terms_multipliable() {
  using lhs_matrix = std::remove_cvref_t<Lhs>;
  using rhs_matrix = std::remove_cvref_t<Rhs>;
  using lhs_row_indexes = typename lhs_matrix::row_indexes;
  using lhs_column_indexes = typename lhs_matrix::column_indexes;
  using rhs_row_indexes = typename rhs_matrix::row_indexes;
  using rhs_column_indexes = typename rhs_matrix::column_indexes;

  bool convertible{true};

  tla::for_constexpr<lhs_matrix::rows>([&convertible](auto i) {
    using lhs_row = tla::product<std::tuple_element_t<i, lhs_row_indexes>,
                                 lhs_column_indexes>;
    tla::for_constexpr<rhs_matrix::columns>([&convertible, i](auto j) {
      using rhs_column =
          tla::product<rhs_row_indexes,
                       std::tuple_element_t<j, rhs_column_indexes>>;
      tla::for_constexpr<lhs_matrix::columns>([&convertible, i, j](auto k) {
        static_cast<void>(i); // Compiler compatibility.
        static_cast<void>(j); // Compiler compatibility.
        convertible &= std::is_convertible_v<
            tla::product<std::tuple_element_t<k, lhs_row>,
                         std::tuple_element_t<k, rhs_column>>,
            tla::product<std::tuple_element_t<0, lhs_row>,
                         std::tuple_element_t<0, rhs_column>>>;
      });
    });
  });

  return convertible;
}

//! @brief Concept of typed matrices whose per-term products can be summed
//! into a matrix product result.
template <typename Lhs, typename Rhs>
concept multipliable_elements = are_terms_multipliable<Lhs, Rhs>();

//! @brief Concept of typed matrices that can be multiplied together.
template <typename Lhs, typename Rhs>
concept multipliable =
    multipliable_shape<Lhs, Rhs> and multipliable_elements<Lhs, Rhs>;

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
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<>;

  return lhs * element{rhs};
}

//! @see operator*(const other auto &, const rank_typed_matrix<0> auto &)
[[nodiscard]] constexpr auto operator*(const rank_typed_matrix<0> auto &lhs,
                                       const other auto &rhs)
  requires scalable_by<
      typename std::remove_cvref_t<decltype(lhs)>::template element<>,
      decltype(rhs)>
{
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using element = typename matrix::template element<>;

  return element{lhs} * rhs;
}

[[nodiscard]] constexpr auto operator*(const rank_typed_matrix<0> auto &lhs,
                                       const rank_typed_matrix<0> auto &rhs) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;
  using lhs_element = typename lhs_matrix::template element<>;
  using rhs_element = typename rhs_matrix::template element<>;

  return lhs_element{lhs} * rhs_element{rhs};
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_PRODUCT_TPP
