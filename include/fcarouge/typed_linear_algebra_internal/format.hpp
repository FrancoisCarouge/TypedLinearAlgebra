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
//!
//! @todo Move the formatter to its header.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes,
          typename Char>
struct std::formatter<fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes>,
                      Char> {
  constexpr auto parse(std::basic_format_parse_context<Char> &parse_context) {
    return parse_context.begin();
  }

  template <typename FormatContext>
  constexpr auto
  format(const fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes> &value,
         FormatContext &format_context) const -> FormatContext::iterator {
    format_context.advance_to(std::format_to(format_context.out(), "["));

    fcarouge::typed_linear_algebra_internal::for_constexpr<
        0, fcarouge::typed_linear_algebra_internal::size<RowIndexes>,
        1>([&value, &format_context](auto i) {
      if (i > 0) {
        format_context.advance_to(std::format_to(format_context.out(), ", "));
      }

      format_context.advance_to(std::format_to(format_context.out(), "["));

      fcarouge::typed_linear_algebra_internal::for_constexpr<
          0, fcarouge::typed_linear_algebra_internal::size<ColumnIndexes>, 1>(
          [&value, &format_context, &i](auto j) {
            if (j > 0) {
              format_context.advance_to(
                  std::format_to(format_context.out(), ", "));
            }

            format_context.advance_to(std::format_to(
                format_context.out(), "{}", value.template at<i, j>()));
          });

      format_context.advance_to(std::format_to(format_context.out(), "]"));
    });

    format_context.advance_to(std::format_to(format_context.out(), "]"));

    return format_context.out();
  }

  template <typename FormatContext>
  constexpr auto
  format(const fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes> &value,
         FormatContext &format_context) const -> FormatContext::iterator
    requires fcarouge::typed_linear_algebra_internal::row<
        fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
  {
    format_context.advance_to(std::format_to(format_context.out(), "["));

    fcarouge::typed_linear_algebra_internal::for_constexpr<
        0, fcarouge::typed_linear_algebra_internal::size<ColumnIndexes>, 1>(
        [&value, &format_context](auto position) {
          if (position > 0) {
            format_context.advance_to(
                std::format_to(format_context.out(), ", "));
          }

          format_context.advance_to(std::format_to(
              format_context.out(), "{}", value.template at<0, position>()));
        });

    format_context.advance_to(std::format_to(format_context.out(), "]"));

    return format_context.out();
  }

  template <typename FormatContext>
  constexpr auto
  format(const fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes> &value,
         FormatContext &format_context) const -> FormatContext::iterator
    requires fcarouge::typed_linear_algebra_internal::singleton<
        fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes>>
  {
    format_context.advance_to(
        std::format_to(format_context.out(), "{}", value.template at<0>()));

    return format_context.out();
  }
};

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_FORMAT_HPP
