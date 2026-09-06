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

#ifndef FCAROUGE_ARMADILLO_HPP
#define FCAROUGE_ARMADILLO_HPP

//! @file
//! @brief Linear algebra facade for the Armadillo third party implementation.
//!
//! @details Supporting compile-time sized matrices and vectors.
//!
//! @note The Armadillo linear algebra is not constexpr-compatible.
//!
//! @note This facade fetches Armadillo for its headers only, in pure
//! header-only mode: no BLAS, LAPACK, ARPACK, or SuperLU is linked. Matrix
//! inversion and `solve()`-backed operations, and therefore matrix division,
//! are unavailable with this backend.

#include "fcarouge/typed_linear_algebra.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <concepts>
#include <cstddef>
#include <format>
#include <tuple>
#include <type_traits>

#include <armadillo>

namespace fcarouge::armadillo {
//! @name Concepts
//! @{

//! @brief An Armadillo algebraic concept.
//!
//! @details Relies on Armadillo's own `is_arma_type` trait: the CRTP base of a
//! nested `Mat<eT>::fixed<R, C>` is `Mat<eT>`, not the `fixed` type itself, so
//! a `std::derived_from<Type, arma::Base<..., Type>>` check would reject every
//! `fixed` matrix.
template <typename Type>
concept is_armadillo = arma::is_arma_type<std::remove_cvref_t<Type>>::value;

//! @brief An Armadillo type with compile-time known dimensions, that is, a
//! `fixed` matrix or vector.
template <typename Type>
concept statically_sized =
    arma::is_Mat_fixed<std::remove_cvref_t<Type>>::value && requires {
      { std::remove_cvref_t<Type>::n_rows } -> std::convertible_to<arma::uword>;
      { std::remove_cvref_t<Type>::n_cols } -> std::convertible_to<arma::uword>;
    };

//! @}

//! @name Types
//! @{

//! @brief Compile-time sized Armadillo matrix.
//!
//! @details Facade for Armadillo implementation compatibility.
//!
//! @tparam Type The matrix element type.
//! @tparam Row The number of rows of the matrix.
//! @tparam Column The number of columns of the matrix.
template <typename Type = double, arma::uword Row = 1, arma::uword Column = 1>
using matrix = typename arma::Mat<Type>::template fixed<Row, Column>;

//! @brief Compile-time sized Armadillo row vector.
template <typename Type = double, arma::uword Column = 1>
using row_vector = typename arma::Mat<Type>::template fixed<1, Column>;

//! @brief Compile-time sized Armadillo column vector.
template <typename Type = double, arma::uword Row = 1>
using column_vector = typename arma::Mat<Type>::template fixed<Row, 1>;

//! @}

} // namespace fcarouge::armadillo

namespace arma {
//! @brief Whole-matrix equality for compile-time sized Armadillo types.
//!
//! @details Armadillo's own `operator==` yields an element-wise `umat`-like
//! expression; this facade mirrors Eigen, whose `operator==` returns a plain
//! `bool`, so the typed matrix equality assertions and the raw backend behave
//! alike. More constrained than Armadillo's operator, so it wins overload
//! resolution for `fixed` operands.
template <typename Lhs, typename Rhs>
  requires fcarouge::armadillo::statically_sized<Lhs> &&
           fcarouge::armadillo::statically_sized<Rhs>
[[nodiscard]] inline auto operator==(const Lhs &lhs, const Rhs &rhs) -> bool {
  static_assert(Lhs::n_rows == Rhs::n_rows && Lhs::n_cols == Rhs::n_cols,
                "Equality requires matrices of the same shape.");

  for (arma::uword k{0}; k < Lhs::n_rows * Lhs::n_cols; ++k) {
    if (lhs[k] != rhs[k]) {
      return false;
    }
  }

  return true;
}

template <typename Lhs, typename Rhs>
  requires fcarouge::armadillo::statically_sized<Lhs> &&
           fcarouge::armadillo::statically_sized<Rhs>
[[nodiscard]] inline auto operator!=(const Lhs &lhs, const Rhs &rhs) -> bool {
  return not(lhs == rhs);
}

//! @brief Get function argument-dependent lookup overload of Armadillo types
//! for structured bindings.
//!
//! @details Row-major linearization, matching the Eigen facade.
//!
//! @todo How should structured bindings be done over a matrix with more than
//! one dimension?
template <std::size_t Index, typename Type>
  requires fcarouge::armadillo::is_armadillo<Type> &&
           fcarouge::armadillo::statically_sized<Type>
[[nodiscard]] constexpr auto get(Type &value) -> typename Type::elem_type & {
  return value.at(Index / Type::n_cols, Index % Type::n_cols);
}

//! @brief Get function argument-dependent lookup overload of Armadillo types
//! for structured bindings.
template <std::size_t Index, typename Type>
  requires fcarouge::armadillo::is_armadillo<Type> &&
           fcarouge::armadillo::statically_sized<Type>
[[nodiscard]] constexpr auto get(const Type &value) ->
    typename Type::elem_type {
  return value.at(Index / Type::n_cols, Index % Type::n_cols);
}
} // namespace arma

//! @brief Opt Armadillo matrices out of the standard range formatter.
//!
//! @details Armadillo matrices are iterable ranges, so without this the
//! library's range formatter and the bracketed matrix formatter below are
//! ambiguous partial specializations of `std::formatter`. Disabling the range
//! interpretation leaves the matrix formatter unambiguously selected. Only
//! relevant, and only available, where the standard library ships the range
//! formatter itself; without it there is no competing specialization and
//! `std::format_kind` does not exist.
//!
//! @note Detected through `_MSC_VER` in addition to the standard
//! `__cpp_lib_format_ranges` feature test macro: MSVC ships an unconditional
//! `std::formatter<Rng, Char>` for any range in `<format>` without yet
//! defining the macro that is supposed to announce it, so relying on the
//! macro alone leaves the opt-out compiled out and the two specializations
//! ambiguous on that standard library.
#if defined(__cpp_lib_format_ranges) || defined(_MSC_VER)
namespace std {
template <typename Type>
  requires fcarouge::armadillo::is_armadillo<Type>
constexpr range_format format_kind<Type> = range_format::disabled;
} // namespace std
#endif

//! @brief Specialization of the standard formatter for the Armadillo matrix.
//!
//! @details Mirrors the bracketed layout of the Eigen facade formatter and of
//! the generic typed matrix formatter: `[[a, b], [c, d]]` for a matrix, `[a,
//! b, c]` for a row vector, and the bare value for a singleton.
//!
//! @note Constrained on the concepts rather than pattern-matched on the alias
//! template, whose `Type`/`Row`/`Column` parameters are not deducible through
//! Armadillo's nested `Mat<eT>::fixed<R, C>` type.
template <typename Type, typename Char>
  requires fcarouge::armadillo::is_armadillo<Type> &&
           fcarouge::armadillo::statically_sized<Type>
struct std::formatter<Type, Char> {
  static constexpr auto
  parse(std::basic_format_parse_context<Char> &parse_context) {
    return parse_context.begin();
  }

  template <typename OutputIterator>
  constexpr auto
  format(const Type &value,
         std::basic_format_context<OutputIterator, Char> &format_context) const
      -> OutputIterator {
    auto output{format_context.out()};

    //! @details Elements are written through `std::to_chars` directly,
    //! rather than through a nested, compile-time-checked
    //! `std::format_to("{}", ...)` call: MSVC's `<format>` immediate function
    //! evaluation for compile-time format string checking has trouble with a
    //! user-defined formatter's own `format()` issuing further checked
    //! format strings. `to_chars`'s default numeric conversion is specified
    //! to produce the same shortest round-trip text as `std::format`, so the
    //! rendered output is unaffected.
    auto write{[&output](const typename Type::elem_type &element) {
      std::array<char, 64> buffer;
      const auto end{
          std::to_chars(buffer.data(), buffer.data() + buffer.size(), element)
              .ptr};
      output = std::copy(buffer.data(), end, output);
    }};

    if constexpr (Type::n_rows == 1 && Type::n_cols == 1) {
      write(value.at(0, 0));
      return output;
    } else if constexpr (Type::n_rows == 1) {
      *output++ = '[';
      for (arma::uword j{0}; j < Type::n_cols; ++j) {
        if (j != 0) {
          *output++ = ',';
          *output++ = ' ';
        }
        write(value.at(0, j));
      }
      *output++ = ']';
      return output;
    } else {
      *output++ = '[';
      for (arma::uword i{0}; i < Type::n_rows; ++i) {
        if (i != 0) {
          *output++ = ',';
          *output++ = ' ';
        }
        *output++ = '[';
        for (arma::uword j{0}; j < Type::n_cols; ++j) {
          if (j != 0) {
            *output++ = ',';
            *output++ = ' ';
          }
          write(value.at(i, j));
        }
        *output++ = ']';
      }
      *output++ = ']';
      return output;
    }
  }
};

//! @brief Tuple size specialization of Armadillo types for structured
//! bindings.
template <typename Type>
  requires fcarouge::armadillo::is_armadillo<Type> &&
           fcarouge::armadillo::statically_sized<Type>
struct std::tuple_size<Type>
    : std::integral_constant<std::size_t, Type::n_rows * Type::n_cols> {};

//! @brief Tuple element specialization of Armadillo types for structured
//! bindings.
template <std::size_t Index, typename Type>
  requires fcarouge::armadillo::is_armadillo<Type> &&
           fcarouge::armadillo::statically_sized<Type> &&
           (Index < std::tuple_size_v<Type>)
struct std::tuple_element<Index, Type> {
  using type = typename Type::elem_type;
};

#endif // FCAROUGE_ARMADILLO_HPP
