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

#ifndef FCAROUGE_LINALG_HPP
#define FCAROUGE_LINALG_HPP

//! @file
//! @brief Indexed-based linear algebra with mp-units with Eigen
//! implementations.

#include "fcarouge/eigen.hpp"
#include "fcarouge/typed_linear_algebra.hpp"
#include "fcarouge/unit.hpp"

namespace fcarouge {
template <typename To, mp_units::Quantity From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr To operator()(const From &value) {
    return value.numerical_value_in(value.unit);
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr To operator()(const From &value) {
    return value * To::reference;
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To &, From &> {
  [[nodiscard]] static constexpr To &operator()(From &value) {
    return reinterpret_cast<To &>(value);
  }
};

template <typename To, mp_units::QuantityPoint From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr To operator()(const From &value) {
    return value.quantity_from_zero().numerical_value_in(value.unit);
  }
};

template <mp_units::QuantityPoint To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr To operator()(const From &value) {
    return {value * To::unit, mp_units::default_point_origin(To::unit)};
  }
};

template <mp_units::QuantityPoint To, typename From>
struct element_caster<To &, From &> {
  [[nodiscard]] static constexpr To &operator()(From &value) {
    return reinterpret_cast<To &>(value);
  }
};

//! @brief Quantity matrix with mp-units and Eigen implementations.
template <typename Representation, typename RowIndexes, typename ColumnIndexes>
using matrix =
    typed_matrix<Eigen::Matrix<Representation, std::tuple_size_v<RowIndexes>,
                               std::tuple_size_v<ColumnIndexes>>,
                 RowIndexes, ColumnIndexes>;

//! @brief Quantity column vector with mp-units and Eigen implementations.
template <typename Representation, typename... Types>
using column_vector =
    typed_column_vector<eigen::column_vector<Representation, sizeof...(Types)>,
                        Types...>;

//! @brief Quantity row vector with mp-units and Eigen implementations.
template <typename Representation, typename... Types>
using row_vector =
    typed_row_vector<eigen::row_vector<Representation, sizeof...(Types)>,
                     Types...>;
} // namespace fcarouge

#endif // FCAROUGE_LINALG_HPP
