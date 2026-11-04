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

#ifndef FCAROUGE_SPARSE_HPP
#define FCAROUGE_SPARSE_HPP

//! @file
//! @brief Compile-time sparsity pattern linear algebra storage.
//!
//! @details The positions of the stored elements are part of the matrix type.
//! Every other element is a structural zero: it occupies no memory and costs no
//! operation. Each operation computes the pattern of its result, and the exact
//! arithmetic terms it needs, at compile time. There is no memory allocation,
//! no exception, no input-output, and every operation is `constexpr`. A pattern
//! without a structural zero behaves as any dense storage.
//!
//! @note A structural zero stays an exact zero: no term exists for it. Unlike
//! a dense computation, a non-finite operand does not propagate through it, for
//! example `0 * infinity` is never computed, and is not a not-a-number value.

#include "fcarouge/typed_linear_algebra.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

#if __STDC_HOSTED__

#include <cmath>

#endif

namespace fcarouge::sparse {
//! @name Patterns
//! @{

//! @brief The index of a position a pattern does not store.
inline constexpr std::size_t absent{std::numeric_limits<std::size_t>::max()};

//! @brief The row and column of a stored element.
struct position {
  std::size_t row{0};
  std::size_t column{0};

  constexpr position() = default;

  constexpr position(std::size_t row_index, std::size_t column_index)
      : row{row_index}, column{column_index} {}

  [[nodiscard]] friend constexpr auto operator==(const position &lhs,
                                                 const position &rhs)
      -> bool = default;

  //! @brief The row-major order.
  [[nodiscard]] friend constexpr auto operator<(const position &lhs,
                                                const position &rhs) -> bool {
    return lhs.row < rhs.row || (lhs.row == rhs.row && lhs.column < rhs.column);
  }
};

//! @brief The positions of the stored elements of a matrix.
//!
//! @details A structural type, usable as a template argument, of strictly
//! ascending row-major positions: two equal patterns are one and the same type
//! of matrix, and a stored element's index is its rank among the positions.
//! Prefer `make_pattern`, `generate`, and the named patterns, which establish
//! that order.
//!
//! @tparam Rows The count of rows of the matrix.
//! @tparam Columns The count of columns of the matrix.
//! @tparam Count The count of stored elements of the matrix.
template <std::size_t Rows, std::size_t Columns, std::size_t Count>
struct pattern {
  std::size_t rows{Rows};
  std::size_t columns{Columns};
  std::size_t count{Count};

  //! @brief The stored positions. An unused position when none is stored:
  //! not every standard library's empty array is a structural type.
  std::array<position, (Count > 0 ? Count : 1)> positions{};

  //! @brief The index of the first stored position not before the given one.
  [[nodiscard]] constexpr auto lower_bound(position value) const
      -> std::size_t {
    std::size_t first{0};
    std::size_t size{Count};

    while (size > 0) {
      const std::size_t half{size / 2};

      if (positions[first + half] < value) {
        first += half + 1;
        size -= half + 1;
      } else {
        size = half;
      }
    }

    return first;
  }

  //! @brief The index of the stored element at the given row and column, or
  //! `absent` for a structural zero.
  [[nodiscard]] constexpr auto index_of(std::size_t row,
                                        std::size_t column) const
      -> std::size_t {
    const position value{row, column};
    const std::size_t index{lower_bound(value)};

    return index < Count && positions[index] == value ? index : absent;
  }

  //! @brief Whether the given row and column is a stored element.
  [[nodiscard]] constexpr auto stores(std::size_t row, std::size_t column) const
      -> bool {
    return index_of(row, column) != absent;
  }
};

//! @brief Whether two patterns store the same positions of the same shape.
template <std::size_t LhsRows, std::size_t LhsColumns, std::size_t LhsCount,
          std::size_t RhsRows, std::size_t RhsColumns, std::size_t RhsCount>
[[nodiscard]] constexpr auto
operator==(const pattern<LhsRows, LhsColumns, LhsCount> &lhs,
           const pattern<RhsRows, RhsColumns, RhsCount> &rhs) -> bool {
  if constexpr (LhsRows != RhsRows || LhsColumns != RhsColumns ||
                LhsCount != RhsCount) {
    return false;
  } else {
    return lhs.positions == rhs.positions;
  }
}

template <typename Type> inline constexpr bool is_pattern{false};

template <std::size_t Rows, std::size_t Columns, std::size_t Count>
inline constexpr bool is_pattern<pattern<Rows, Columns, Count>>{true};

//! @brief Concept of a compile-time sparsity pattern.
template <typename Type>
concept sparsity_pattern = is_pattern<std::remove_cv_t<Type>>;

//! @brief Called by a pattern's invalid construction: never defined, its call
//! ends the constant evaluation, and is the name the compiler's diagnostic
//! reports.
void position_out_of_range_or_duplicate();

//! @brief The pattern of an empty matrix.
template <std::size_t Rows, std::size_t Columns>
[[nodiscard]] consteval auto make_pattern() -> pattern<Rows, Columns, 0> {
  return {};
}

//! @brief The pattern of the given positions, in any order.
//!
//! @details A position out of the matrix, or a duplicate position, does not
//! compile.
template <std::size_t Rows, std::size_t Columns, std::size_t Count>
[[nodiscard]] consteval auto make_pattern(const position (&positions)[Count])
    -> pattern<Rows, Columns, Count> {
  pattern<Rows, Columns, Count> result;

  for (std::size_t i{0}; i < Count; ++i) {
    const position value{positions[i]};
    std::size_t j{i};

    if (value.row >= Rows || value.column >= Columns) {
      position_out_of_range_or_duplicate();
    }

    for (; j > 0 && value < result.positions[j - 1]; --j) {
      result.positions[j] = result.positions[j - 1];
    }

    if (j > 0 && value == result.positions[j - 1]) {
      position_out_of_range_or_duplicate();
    }

    result.positions[j] = value;
  }

  return result;
}

//! @brief The pattern of the positions satisfying a predicate.
//!
//! @tparam Predicate A stateless predicate of a row and a column.
template <std::size_t Rows, std::size_t Columns, auto Predicate>
[[nodiscard]] consteval auto generate() {
  constexpr std::size_t count{[] {
    std::size_t result{0};

    for (std::size_t row{0}; row < Rows; ++row) {
      for (std::size_t column{0}; column < Columns; ++column) {
        result += Predicate(row, column) ? std::size_t{1} : std::size_t{0};
      }
    }

    return result;
  }()};

  pattern<Rows, Columns, count> result;

  for (std::size_t row{0}, index{0}; row < Rows; ++row) {
    for (std::size_t column{0}; column < Columns; ++column) {
      if (Predicate(row, column)) {
        result.positions[index] = {row, column};
        ++index;
      }
    }
  }

  return result;
}

//! @brief Every position is stored: a matrix without structural zeros.
template <std::size_t Rows, std::size_t Columns>
inline constexpr auto dense{
    generate<Rows, Columns,
             [](std::size_t, std::size_t) -> bool { return true; }>()};

//! @brief No position is stored: the zero matrix.
template <std::size_t Rows, std::size_t Columns>
inline constexpr auto zero{make_pattern<Rows, Columns>()};

//! @brief The positions of the `Lower` sub-diagonals, the main diagonal, and
//! the `Upper` super-diagonals.
template <std::size_t Rows, std::size_t Columns, std::size_t Lower,
          std::size_t Upper>
inline constexpr auto banded{
    generate<Rows, Columns, [](std::size_t row, std::size_t column) -> bool {
      return row <= column + Lower && column <= row + Upper;
    }>()};

//! @brief The main diagonal of a square matrix.
template <std::size_t Size>
inline constexpr auto diagonal{banded<Size, Size, 0, 0>};

//! @brief The three central diagonals of a square matrix.
template <std::size_t Size>
inline constexpr auto tridiagonal{banded<Size, Size, 1, 1>};

//! @brief The main diagonal and below of a square matrix.
template <std::size_t Size>
inline constexpr auto lower_triangular{banded<Size, Size, Size, 0>};

//! @brief The main diagonal and above of a square matrix.
template <std::size_t Size>
inline constexpr auto upper_triangular{banded<Size, Size, 0, Size>};

//! @brief A growing list of positions, of a known maximum size, from which the
//! pattern of an operation's result is sized.
template <std::size_t Capacity> struct position_list {
  std::array<position, (Capacity > 0 ? Capacity : 1)> items{};
  std::size_t size{0};

  constexpr void push(position value) {
    items[size] = value;
    ++size;
  }
};

template <std::size_t Rows, std::size_t Columns, auto List>
[[nodiscard]] consteval auto to_pattern() -> pattern<Rows, Columns, List.size> {
  pattern<Rows, Columns, List.size> result;

  for (std::size_t i{0}; i < List.size; ++i) {
    result.positions[i] = List.items[i];
  }

  return result;
}

//! @brief Merge two patterns, keeping the positions of both, or only those in
//! common, in ascending order.
template <auto Lhs, auto Rhs, bool KeepAll>
[[nodiscard]] consteval auto merge_positions() {
  position_list<Lhs.count + Rhs.count> list;
  std::size_t i{0};
  std::size_t j{0};

  while (i < Lhs.count || j < Rhs.count) {
    if (j == Rhs.count ||
        (i < Lhs.count && Lhs.positions[i] < Rhs.positions[j])) {
      if constexpr (KeepAll) {
        list.push(Lhs.positions[i]);
      }
      ++i;
    } else if (i == Lhs.count || Rhs.positions[j] < Lhs.positions[i]) {
      if constexpr (KeepAll) {
        list.push(Rhs.positions[j]);
      }
      ++j;
    } else {
      list.push(Lhs.positions[i]);
      ++i;
      ++j;
    }
  }

  return list;
}

//! @brief The positions stored by either pattern of the same shape: the pattern
//! of a sum or a difference.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
inline constexpr auto pattern_union{
    to_pattern<Lhs.rows, Lhs.columns, merge_positions<Lhs, Rhs, true>()>()};

//! @brief The positions stored by both patterns of the same shape.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
inline constexpr auto pattern_intersection{
    to_pattern<Lhs.rows, Lhs.columns, merge_positions<Lhs, Rhs, false>()>()};

//! @brief The transposed pattern.
//!
//! @details The positions are bucketed by column: ascending rows within each
//! column are already what the transposed row-major order requires.
template <sparsity_pattern auto Pattern>
inline constexpr auto pattern_transpose{[] {
  pattern<Pattern.columns, Pattern.rows, Pattern.count> result;
  std::size_t index{0};

  for (std::size_t column{0}; column < Pattern.columns; ++column) {
    for (std::size_t i{0}; i < Pattern.count; ++i) {
      if (Pattern.positions[i].column == column) {
        result.positions[index] = {column, Pattern.positions[i].row};
        ++index;
      }
    }
  }

  return result;
}()};

//! @brief The count of elementary products of a matrix product.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
inline constexpr std::size_t product_term_count{[] {
  std::size_t result{0};

  for (std::size_t i{0}; i < Lhs.count; ++i) {
    const std::size_t row{Lhs.positions[i].column};

    result += Rhs.lower_bound({row + 1, 0}) - Rhs.lower_bound({row, 0});
  }

  return result;
}()};

//! @brief The positions reached by at least one term of a matrix product.
template <auto Lhs, auto Rhs> [[nodiscard]] consteval auto product_positions() {
  position_list<product_term_count<Lhs, Rhs>> list;

  for (std::size_t row{0}; row < Lhs.rows; ++row) {
    std::array<bool, Rhs.columns> reached{};
    const std::size_t end{Lhs.lower_bound({row + 1, 0})};

    for (std::size_t i{Lhs.lower_bound({row, 0})}; i < end; ++i) {
      const std::size_t inner{Lhs.positions[i].column};
      const std::size_t inner_end{Rhs.lower_bound({inner + 1, 0})};

      for (std::size_t j{Rhs.lower_bound({inner, 0})}; j < inner_end; ++j) {
        reached[Rhs.positions[j].column] = true;
      }
    }

    for (std::size_t column{0}; column < Rhs.columns; ++column) {
      if (reached[column]) {
        list.push({row, column});
      }
    }
  }

  return list;
}

//! @brief The structural pattern of a matrix product.
//!
//! @details A position is stored whenever a term reaches it, even if the values
//! of the terms may cancel each other out at run time.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.columns == Rhs.rows)
inline constexpr auto pattern_product{
    to_pattern<Lhs.rows, Rhs.columns, product_positions<Lhs, Rhs>()>()};

//! @brief Whether a pattern stores every position of another.
template <sparsity_pattern auto Inner, sparsity_pattern auto Outer>
inline constexpr bool is_subset_of{[] {
  if constexpr (Inner.rows != Outer.rows or Inner.columns != Outer.columns) {
    return false;
  } else {
    return pattern_intersection<Inner, Outer>.count == Inner.count;
  }
}()};

//! @}

//! @name Operation Tables
//! @{

//! @brief For each position of `To`, its stored element index in `From`, or
//! `absent`.
template <auto From, auto To>
inline constexpr auto projection{[] {
  std::array<std::size_t, To.count> result{};

  for (std::size_t i{0}; i < To.count; ++i) {
    result[i] = From.index_of(To.positions[i].row, To.positions[i].column);
  }

  return result;
}()};

//! @brief An elementary product of a matrix product, or of a dot product.
struct term {
  std::size_t result{0};
  std::size_t lhs{0};
  std::size_t rhs{0};
};

//! @brief The elementary products of a matrix product.
template <auto Lhs, auto Rhs>
inline constexpr auto product_terms{[] {
  constexpr auto result_pattern{pattern_product<Lhs, Rhs>};
  std::array<term, product_term_count<Lhs, Rhs>> result{};
  std::size_t index{0};

  for (std::size_t i{0}; i < Lhs.count; ++i) {
    const auto [row, inner]{Lhs.positions[i]};
    const std::size_t inner_end{Rhs.lower_bound({inner + 1, 0})};

    for (std::size_t j{Rhs.lower_bound({inner, 0})}; j < inner_end; ++j) {
      result[index] = {result_pattern.index_of(row, Rhs.positions[j].column), i,
                       j};
      ++index;
    }
  }

  return result;
}()};

//! @brief The position along a vector of the given stored element.
template <auto Pattern>
[[nodiscard]] constexpr auto vector_index(std::size_t i) -> std::size_t {
  return (Pattern.positions[i].row * Pattern.columns) +
         Pattern.positions[i].column;
}

//! @brief The elementary products of the dot product of two vectors, of the
//! same size and of any orientation.
template <auto Lhs, auto Rhs>
inline constexpr auto dot_terms{[] {
  constexpr std::size_t count{[] {
    std::size_t result{0};

    for (std::size_t i{0}, j{0}; i < Lhs.count && j < Rhs.count;) {
      if (vector_index<Lhs>(i) < vector_index<Rhs>(j)) {
        ++i;
      } else if (vector_index<Rhs>(j) < vector_index<Lhs>(i)) {
        ++j;
      } else {
        ++result;
        ++i;
        ++j;
      }
    }

    return result;
  }()};
  std::array<term, count> result{};
  std::size_t index{0};

  for (std::size_t i{0}, j{0}; i < Lhs.count && j < Rhs.count;) {
    if (vector_index<Lhs>(i) < vector_index<Rhs>(j)) {
      ++i;
    } else if (vector_index<Rhs>(j) < vector_index<Lhs>(i)) {
      ++j;
    } else {
      result[index] = {index, i, j};
      ++index;
      ++i;
      ++j;
    }
  }

  return result;
}()};

//! @}

//! @name Types
//! @{

//! @brief The square root, without a hosted implementation or at compile time
//! by Newton's method, descending from above the root.
template <std::floating_point Type>
[[nodiscard]] constexpr auto square_root(Type value) -> Type {
#if __STDC_HOSTED__
  if !consteval {
    return std::sqrt(value);
  }
#endif

  // A zero, a negative, a not-a-number, or an infinite value has no iteration.
  if (not(value > Type{0}) || value > std::numeric_limits<Type>::max()) {
    return value;
  }

  Type estimate{value < Type{1} ? Type{1} : value};

  for (Type next{(estimate + (value / estimate)) / Type{2}}; next < estimate;
       next = (estimate + (value / estimate)) / Type{2}) {
    estimate = next;
  }

  return estimate;
}

//! @brief The sums of the rows or the columns of a matrix.
//!
//! @details Eagerly computed. The `sum` and `maxCoeff` members match Eigen's
//! naming, through which the norms reach this backend's storage.
template <typename Type, std::size_t Size> struct line_sums {
  std::array<Type, Size> sums{};

  [[nodiscard]] constexpr auto sum() const -> line_sums { return *this; }

  [[nodiscard]] constexpr auto maxCoeff() const -> Type {
    Type result{sums[0]};

    for (const Type &value : sums) {
      result = result < value ? value : result;
    }

    return result;
  }
};

//! @brief Compile-time sparsity pattern matrix storage.
//!
//! @details Stores only the elements of the pattern, in row-major order, in
//! place. A read of a structural zero yields a zero value. A pattern with a
//! structural zero only has an element reference for a compile-time stored
//! position, so no write can reach a structural zero and no value can be lost:
//! its stored elements are written with the typed matrix's `at<i, j>(value)` or
//! `sparse::make_matrix`.
//!
//! @tparam Type The matrix element type.
//! @tparam Pattern The sparsity pattern of the matrix.
template <typename Type, sparsity_pattern auto Pattern> struct matrix {
  //! @brief The sparsity pattern of the matrix.
  static constexpr auto structure{Pattern};

  //! @brief The stored elements, zero-initialized, in row-major order.
  std::array<Type, Pattern.count> values{};

  //! @brief Construct a zero matrix.
  constexpr matrix() = default;

  //! @brief Construct from every stored element, in row-major order.
  template <std::same_as<Type>... Values>
    requires(sizeof...(Values) == Pattern.count and Pattern.count > 0)
  constexpr explicit matrix(const Values &...stored) : values{stored...} {}

  //! @brief Convert from a matrix whose stored positions are all stored by this
  //! pattern, which cannot lose an element.
  template <sparsity_pattern auto Other>
    requires is_subset_of<Other, Pattern>
  constexpr matrix(const matrix<Type, Other> &other) {
    for (std::size_t i{0}; i < Pattern.count; ++i) {
      const std::size_t index{projection<Other, Pattern>[i]};

      if (index != absent) {
        values[i] = other.values[index];
      }
    }
  }

  //! @brief Read the element at the given position, by value.
  [[nodiscard]] constexpr auto operator[](std::size_t row,
                                          std::size_t column) const -> Type {
    const std::size_t index{Pattern.index_of(row, column)};

    return index == absent ? Type{} : values[index];
  }

  //! @brief The element at the given position, of a pattern without a
  //! structural zero: every position is addressable, as in any dense storage.
  [[nodiscard]] constexpr auto operator[](std::size_t row, std::size_t column)
      -> Type &
    requires(Pattern.count == Pattern.rows * Pattern.columns)
  {
    return values[(row * Pattern.columns) + column];
  }

  //! @brief The stored element at the given compile-time position. A structural
  //! zero has nothing to refer to, and no such access.
  template <std::size_t Row, std::size_t Column>
    requires(Pattern.stores(Row, Column))
  [[nodiscard]] constexpr auto
  operator[](std::integral_constant<std::size_t, Row> /*row*/,
             std::integral_constant<std::size_t, Column> /*column*/) -> Type & {
    static constexpr std::size_t index{Pattern.index_of(Row, Column)};

    return values[index];
  }

  //! @brief The transposed matrix.
  [[nodiscard]] constexpr auto transpose() const {
    constexpr auto transposed{pattern_transpose<Pattern>};
    constexpr auto indexes{[&transposed] {
      std::array<std::size_t, Pattern.count> result{};

      for (std::size_t i{0}; i < Pattern.count; ++i) {
        result[i] = Pattern.index_of(transposed.positions[i].column,
                                     transposed.positions[i].row);
      }

      return result;
    }()};
    matrix<Type, transposed> result;

    for (std::size_t i{0}; i < Pattern.count; ++i) {
      result.values[i] = values[indexes[i]];
    }

    return result;
  }

  //! @brief The sum of the products of the elements stored by both vectors.
  template <sparsity_pattern auto Other>
    requires((Pattern.rows == 1 or Pattern.columns == 1) and
             (Other.rows == 1 or Other.columns == 1) and
             Pattern.rows * Pattern.columns == Other.rows * Other.columns)
  [[nodiscard]] constexpr auto dot(const matrix<Type, Other> &other) const
      -> Type {
    Type result{};

    for (const term &element : dot_terms<Pattern, Other>) {
      result += values[element.lhs] * other.values[element.rhs];
    }

    return result;
  }

  //! @brief The Frobenius norm.
  [[nodiscard]] constexpr auto norm() const -> Type
    requires std::floating_point<Type>
  {
    Type result{};

    for (const Type &value : values) {
      result += value * value;
    }

    return square_root(result);
  }

  //! @brief The element-wise absolute values.
  [[nodiscard]] constexpr auto cwiseAbs() const -> matrix {
    matrix result{*this};

    for (Type &value : result.values) {
      value = value < Type{} ? -value : value;
    }

    return result;
  }

  //! @brief The sums of the columns.
  [[nodiscard]] constexpr auto colwise() const
      -> line_sums<Type, Pattern.columns> {
    line_sums<Type, Pattern.columns> result;

    for (std::size_t i{0}; i < Pattern.count; ++i) {
      result.sums[Pattern.positions[i].column] += values[i];
    }

    return result;
  }

  //! @brief The sums of the rows.
  [[nodiscard]] constexpr auto rowwise() const
      -> line_sums<Type, Pattern.rows> {
    line_sums<Type, Pattern.rows> result;

    for (std::size_t i{0}; i < Pattern.count; ++i) {
      result.sums[Pattern.positions[i].row] += values[i];
    }

    return result;
  }
};

//! @}

//! @name Operators
//! @{

//! @brief Combine the elements stored by either matrix of the same shape.
//!
//! @param both The operation on the elements stored by both matrices.
//! @param rhs_only The operation on the elements stored by the rhs matrix only.
//! An element stored by the lhs matrix only is kept as-is.
template <typename Type, auto Lhs, auto Rhs>
[[nodiscard]] constexpr auto combine(const matrix<Type, Lhs> &lhs,
                                     const matrix<Type, Rhs> &rhs,
                                     const auto &both, const auto &rhs_only) {
  constexpr auto combined{pattern_union<Lhs, Rhs>};
  matrix<Type, combined> result;

  for (std::size_t i{0}; i < combined.count; ++i) {
    const std::size_t lhs_index{projection<Lhs, combined>[i]};
    const std::size_t rhs_index{projection<Rhs, combined>[i]};

    if (lhs_index == absent) {
      result.values[i] = rhs_only(rhs.values[rhs_index]);
    } else if (rhs_index == absent) {
      result.values[i] = lhs.values[lhs_index];
    } else {
      result.values[i] = both(lhs.values[lhs_index], rhs.values[rhs_index]);
    }
  }

  return result;
}

//! @brief Addition, whose pattern is the union of the terms' patterns.
template <typename Type, auto Lhs, auto Rhs>
[[nodiscard]] constexpr auto operator+(const matrix<Type, Lhs> &lhs,
                                       const matrix<Type, Rhs> &rhs) {
  return combine(
      lhs, rhs,
      [](const Type &left, const Type &right) -> Type { return left + right; },
      [](const Type &right) -> Type { return right; });
}

//! @brief Subtraction, whose pattern is the union of the terms' patterns.
template <typename Type, auto Lhs, auto Rhs>
[[nodiscard]] constexpr auto operator-(const matrix<Type, Lhs> &lhs,
                                       const matrix<Type, Rhs> &rhs) {
  return combine(
      lhs, rhs,
      [](const Type &left, const Type &right) -> Type { return left - right; },
      [](const Type &right) -> Type { return -right; });
}

//! @brief Negation, of the same pattern.
template <typename Type, auto Pattern>
[[nodiscard]] constexpr auto operator-(const matrix<Type, Pattern> &value)
    -> matrix<Type, Pattern> {
  matrix<Type, Pattern> result{value};

  for (Type &element : result.values) {
    element = -element;
  }

  return result;
}

//! @brief Matrix multiplication, whose pattern is the structural product.
template <typename Type, auto Lhs, auto Rhs>
  requires(Lhs.columns == Rhs.rows)
[[nodiscard]] constexpr auto operator*(const matrix<Type, Lhs> &lhs,
                                       const matrix<Type, Rhs> &rhs) {
  matrix<Type, pattern_product<Lhs, Rhs>> result;

  for (const term &element : product_terms<Lhs, Rhs>) {
    result.values[element.result] +=
        lhs.values[element.lhs] * rhs.values[element.rhs];
  }

  return result;
}

//! @brief Scaling by a scalar factor, on the right, of the same pattern.
template <typename Type, auto Pattern>
[[nodiscard]] constexpr auto operator*(const matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> matrix<Type, Pattern> {
  matrix<Type, Pattern> result{lhs};

  for (Type &element : result.values) {
    element = element * rhs;
  }

  return result;
}

//! @brief Scaling by a scalar factor, on the left, of the same pattern.
template <typename Type, auto Pattern>
[[nodiscard]] constexpr auto operator*(const std::type_identity_t<Type> &lhs,
                                       const matrix<Type, Pattern> &rhs)
    -> matrix<Type, Pattern> {
  matrix<Type, Pattern> result{rhs};

  for (Type &element : result.values) {
    element = lhs * element;
  }

  return result;
}

//! @brief Scaling by the inverse of a scalar divisor, of the same pattern.
//!
//! @note A zero divisor leaves every structural zero an exact zero, see the
//! file's note on non-finite values.
template <typename Type, auto Pattern>
[[nodiscard]] constexpr auto operator/(const matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> matrix<Type, Pattern> {
  matrix<Type, Pattern> result{lhs};

  for (Type &element : result.values) {
    element = element / rhs;
  }

  return result;
}

//! @brief Element-wise equality of two matrices of the same shape, of any
//! pattern: a structural zero equals a stored zero.
template <typename Type, auto Lhs, auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator==(const matrix<Type, Lhs> &lhs,
                                        const matrix<Type, Rhs> &rhs) -> bool {
  constexpr auto combined{pattern_union<Lhs, Rhs>};

  for (std::size_t i{0}; i < combined.count; ++i) {
    const auto [row, column]{combined.positions[i]};

    if (not(lhs[row, column] == rhs[row, column])) {
      return false;
    }
  }

  return true;
}

//! @}
} // namespace fcarouge::sparse

namespace fcarouge {
//! @name Types
//! @{

//! @brief Typed matrix of a compile-time sparsity pattern.
//!
//! @tparam Representation The stored element representation type.
//! @tparam Pattern The sparsity pattern, of the indexes' shape.
//! @tparam RowIndexes The row indexes, by default of the representation type.
//! @tparam ColumnIndexes The column indexes, by default of the representation
//! type.
template <typename Representation, sparse::sparsity_pattern auto Pattern,
          typename RowIndexes = tla::tuple_n_type<Representation, Pattern.rows>,
          typename ColumnIndexes =
              tla::tuple_n_type<Representation, Pattern.columns>>
  requires(Pattern.rows == std::tuple_size_v<RowIndexes> and
           Pattern.columns == std::tuple_size_v<ColumnIndexes>)
using sparse_matrix = typed_matrix<sparse::matrix<Representation, Pattern>,
                                   RowIndexes, ColumnIndexes>;

//! @brief Typed column vector of a compile-time sparsity pattern.
template <typename Representation, sparse::sparsity_pattern auto Pattern,
          typename... Types>
  requires(Pattern.rows == sizeof...(Types) and Pattern.columns == 1)
using sparse_column_vector =
    typed_column_vector<sparse::matrix<Representation, Pattern>, Types...>;

//! @brief Typed row vector of a compile-time sparsity pattern.
template <typename Representation, sparse::sparsity_pattern auto Pattern,
          typename... Types>
  requires(Pattern.rows == 1 and Pattern.columns == sizeof...(Types))
using sparse_row_vector =
    typed_row_vector<sparse::matrix<Representation, Pattern>, Types...>;

//! @}
} // namespace fcarouge

namespace fcarouge::sparse {
//! @name Element Construction
//! @{

//! @brief The value of a structural zero in `make_matrix`.
struct structural_zero_t {};

//! @brief The value of a structural zero in `make_matrix`.
inline constexpr structural_zero_t structural_zero;

//! @brief Whether a value may be given at a linear position of a typed matrix:
//! a structural zero only takes `structural_zero`, a stored element only takes
//! an implicitly convertible value.
//!
//! @details Partially specialized variable templates rather than immediately
//! invoked templated lambdas, for the constraints' evaluation by MSVC.
template <typename Value, typename TypedMatrix, std::size_t Position,
          bool Stored = TypedMatrix::matrix::structure.stores(
              Position / TypedMatrix::columns, Position % TypedMatrix::columns)>
inline constexpr bool is_given_at{std::is_convertible_v<
    Value, tla::element_at<TypedMatrix, Position / TypedMatrix::columns,
                           Position % TypedMatrix::columns>>};

template <typename Value, typename TypedMatrix, std::size_t Position>
inline constexpr bool is_given_at<Value, TypedMatrix, Position, false>{
    std::same_as<Value, structural_zero_t>};

template <typename TypedMatrix, typename Positions, typename... Values>
inline constexpr bool are_given_at{false};

template <typename TypedMatrix, std::size_t... Positions, typename... Values>
inline constexpr bool
    are_given_at<TypedMatrix, std::index_sequence<Positions...>, Values...>{
        (is_given_at<Values, TypedMatrix, Positions> and ...)};

template <typename TypedMatrix, typename... Values>
inline constexpr bool are_given{
    are_given_at<TypedMatrix, std::index_sequence_for<Values...>, Values...>};

//! @brief A typed matrix from a value at every position, in row-major order.
//!
//! @details Each position is known at compile time, and the type of each value
//! says what it is: a value for a structural zero, even a zero, or
//! `structural_zero` for a stored element, does not compile. A run-time value
//! needs no run-time check.
//!
//! @tparam TypedMatrix The typed matrix over a sparse storage to construct.
template <typename TypedMatrix, typename... Values>
  requires(sizeof...(Values) == TypedMatrix::rows * TypedMatrix::columns) and
          are_given<TypedMatrix, Values...>
[[nodiscard]] constexpr auto make_matrix(const Values &...values)
    -> TypedMatrix {
  TypedMatrix result;

  [&result, &values...]<std::size_t... Positions>(
      std::index_sequence<Positions...>) -> void {
    const auto write{[&result]<std::size_t Position>(
                         const auto &value) -> void {
      if constexpr (not std::same_as<std::remove_cvref_t<decltype(value)>,
                                     structural_zero_t>) {
        static constexpr std::size_t row{Position / TypedMatrix::columns};
        static constexpr std::size_t column{Position % TypedMatrix::columns};
        static constexpr std::size_t index{
            TypedMatrix::matrix::structure.index_of(row, column)};

        using element = tla::element_at<TypedMatrix, row, column>;

        // Not braces: a conversion to the element type may narrow.
        const element converted(value);

        result.data().values[index] =
            cast<typename TypedMatrix::underlying, element>(converted);
      }
    }};

    (write.template operator()<Positions>(values), ...);
  }(std::index_sequence_for<Values...>{});

  return result;
}

//! @}
} // namespace fcarouge::sparse

#endif // FCAROUGE_SPARSE_HPP
