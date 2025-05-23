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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_FORMAT_HPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_FORMAT_HPP

//! @file
//! @brief Formatting support for the typed matrix.

#include "utility.hpp"

#include <cstddef>
#include <format>

//! @brief Specialization of the standard formatter for the typed matrix.
template <fcarouge::typed_linear_algebra_internal::typed_matrix TypedMatrix,
          typename Char>
struct std::formatter<TypedMatrix, Char> {
  constexpr auto parse(std::basic_format_parse_context<Char> &parse_context) {
    return parse_context.begin();
  }

  template <typename OutputIterator>
  constexpr auto
  format(const TypedMatrix &value,
         std::basic_format_context<OutputIterator, Char> &format_context) const
      -> OutputIterator {
    format_context.advance_to(std::format_to(format_context.out(), "["));

    for (std::size_t i{0}; i < TypedMatrix::rows; ++i) {
      if (i > 0) {
        format_context.advance_to(std::format_to(format_context.out(), ", "));
      }

      format_context.advance_to(std::format_to(format_context.out(), "["));

      for (std::size_t j{0}; j < TypedMatrix::columns; ++j) {
        if (j > 0) {
          format_context.advance_to(std::format_to(format_context.out(), ", "));
        }

        format_context.advance_to(
            std::format_to(format_context.out(), "{}", value.data(i, j)));
      }

      format_context.advance_to(std::format_to(format_context.out(), "]"));
    }

    format_context.advance_to(std::format_to(format_context.out(), "]"));

    return format_context.out();
  }

  template <typename OutputIterator>
  constexpr auto
  format(const TypedMatrix &value,
         std::basic_format_context<OutputIterator, Char> &format_context) const
      -> OutputIterator
    requires fcarouge::typed_linear_algebra_internal::row<TypedMatrix>
  {
    format_context.advance_to(std::format_to(format_context.out(), "["));

    for (std::size_t j{0}; j < TypedMatrix::columns; ++j) {
      if (j > 0) {
        format_context.advance_to(std::format_to(format_context.out(), ", "));
      }

      format_context.advance_to(
          std::format_to(format_context.out(), "{}", value.data(0, j)));
    }

    format_context.advance_to(std::format_to(format_context.out(), "]"));

    return format_context.out();
  }

  template <typename OutputIterator>
  constexpr auto
  format(const TypedMatrix &value,
         std::basic_format_context<OutputIterator, Char> &format_context) const
      -> OutputIterator
    requires fcarouge::typed_linear_algebra_internal::singleton<TypedMatrix>
  {
    format_context.advance_to(
        std::format_to(format_context.out(), "{}", value.data(0, 0)));

    return format_context.out();
  }
};

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_FORMAT_HPP
