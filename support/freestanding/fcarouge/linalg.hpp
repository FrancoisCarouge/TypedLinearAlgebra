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

#ifndef FCAROUGE_LINALG_HPP
#define FCAROUGE_LINALG_HPP

//! @file
//! @brief Scalar type linear algebra with a freestanding implementation.
//!
//! @details A fixed-size, row-major matrix storage over a standard array,
//! relying only on the freestanding subset of the standard library: no memory
//! allocation, no exception, no input-output. None of the third party linear
//! algebra backends supports a freestanding implementation.

#include "fcarouge/typed_linear_algebra.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <utility>

namespace fcarouge::freestanding {
//! @name Types
//! @{

//! @brief Compile-time sized matrix storage.
//!
//! @details Supporting the element access, the addition, the subtraction, the
//! negation, the matrix multiplication, and the scaling by a scalar factor.
//!
//! @tparam Type The matrix element type.
//! @tparam Row The number of rows of the matrix.
//! @tparam Column The number of columns of the matrix.
template <typename Type = double, std::size_t Row = 1, std::size_t Column = 1>
struct matrix {
  //! @brief The row-major, zero-initialized elements.
  std::array<Type, Row * Column> elements{};

  //! @brief Access the element at the given row and column.
  [[nodiscard]] constexpr auto operator[](this auto &&self, std::size_t row,
                                          std::size_t column) -> auto & {
    return self.elements[(row * Column) + column];
  }

  //! @brief Element-wise addition of two matrices of the same shape.
  [[nodiscard]] friend constexpr auto operator+(const matrix &lhs,
                                                const matrix &rhs) -> matrix {
    return transform(std::plus{}, lhs, rhs);
  }

  //! @brief Element-wise subtraction of two matrices of the same shape.
  [[nodiscard]] friend constexpr auto operator-(const matrix &lhs,
                                                const matrix &rhs) -> matrix {
    return transform(std::minus{}, lhs, rhs);
  }

  //! @brief Element-wise negation of a matrix.
  [[nodiscard]] friend constexpr auto operator-(const matrix &value) -> matrix {
    return transform(std::negate{}, value);
  }

  //! @brief Scaling of a matrix by a scalar factor, on the right.
  [[nodiscard]] friend constexpr auto operator*(const matrix &lhs,
                                                const Type &rhs) -> matrix {
    return transform(
        [&rhs](const Type &element) -> Type { return element * rhs; }, lhs);
  }

  //! @brief Scaling of a matrix by a scalar factor, on the left.
  [[nodiscard]] friend constexpr auto operator*(const Type &lhs,
                                                const matrix &rhs) -> matrix {
    return transform(
        [&lhs](const Type &element) -> Type { return lhs * element; }, rhs);
  }

  //! @brief Matrix multiplication.
  template <std::size_t Other>
  [[nodiscard]] friend constexpr auto
  operator*(const matrix &lhs, const matrix<Type, Column, Other> &rhs)
      -> matrix<Type, Row, Other> {
    matrix<Type, Row, Other> result;

    for (std::size_t i{0}; i < Row; ++i) {
      for (std::size_t j{0}; j < Other; ++j) {
        for (std::size_t k{0}; k < Column; ++k) {
          result[i, j] += lhs[i, k] * rhs[k, j];
        }
      }
    }

    return result;
  }

private:
  //! @brief The function applied to the elements of the given matrices at the
  //! given position.
  template <std::size_t Position>
  [[nodiscard]] static constexpr auto transform_at(const auto &function,
                                                   const auto &...values)
      -> Type {
    return function(std::get<Position>(values.elements)...);
  }

  //! @brief The matrix of the function applied to the elements of the given
  //! matrices at each position.
  [[nodiscard]] static constexpr auto transform(const auto &function,
                                                const auto &...values)
      -> matrix {
    return [&]<std::size_t... Positions>(
               std::index_sequence<Positions...>) -> matrix {
      return {{transform_at<Positions>(function, values...)...}};
    }(std::make_index_sequence<Row * Column>{});
  }
};

//! @}
} // namespace fcarouge::freestanding

namespace fcarouge {
//! @name Types
//! @{

//! @brief Scalar type matrix with a freestanding implementation.
template <typename Type = double, std::size_t Row = 1, std::size_t Column = 1>
using matrix =
    typed_matrix<freestanding::matrix<Type, Row, Column>,
                 tla::tuple_n_type<Type, Row>, tla::tuple_n_type<Type, Column>>;

//! @brief Scalar type column vector with a freestanding implementation.
template <typename Type = double, std::size_t Row = 1>
using column_vector = matrix<Type, Row, 1>;

//! @brief Scalar type row vector with a freestanding implementation.
template <typename Type = double, std::size_t Column = 1>
using row_vector = matrix<Type, 1, Column>;

//! @}
} // namespace fcarouge

#endif // FCAROUGE_LINALG_HPP
