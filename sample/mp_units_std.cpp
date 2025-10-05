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

//! @file
//! @brief Unit safe linear algebra with mp-units and Eigen.
//!
//! @details Demonstrate a variety of linear algebra operations with mp-units
//! and Eigen. This library composes Eigen as the matrix' linear algebra backend
//! with index typed as mp-units types. This sample explicitly uses double
//! precision floating point numbers. This sample uses Eigen linear algebra as
//! the linear algebra backend. This sample uses mp-units types for the strongly
//! typed units.

// #include "fcarouge/linalg.hpp"
#include "fcarouge/typed_linear_algebra.hpp"

#include <format>
#include <print>
#include <tuple>
#include <type_traits>

#include <mp-units/framework/quantity.h>
#include <mp-units/framework/quantity_point.h>
#include <mp-units/math.h>
#include <mp-units/systems/isq/thermodynamics.h>
#include <mp-units/systems/si.h>

namespace fcarouge {
// Teach the typed linear algebra library how to convert Eigen' underlying
// scalar types to and from mp-units' types.
template <typename To, mp_units::Quantity From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(From value) const {
    return value.numerical_value_in(value.unit);
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(From value) const {
    return value * To::reference;
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To &, From &> {
  [[nodiscard]] constexpr To &operator()(From &value) const {
    return reinterpret_cast<To &>(value);
  }
};

template <typename To, mp_units::Reference From>
struct element_caster<To, From> {
  [[nodiscard]] constexpr To operator()(From) const { return 1.; }
};

namespace sample {
namespace {
// Set up heterogenously unit typed linear algebra types.
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

template <mp_units::Reference auto QuantityReference>
using quantity_point =
    mp_units::quantity_point<QuantityReference,
                             mp_units::default_point_origin(QuantityReference),
                             representation>;

// template <typename RowIndexes, typename ColumnIndexes>
// using matrix =
//     typed_matrix<Eigen::Matrix<representation, std::tuple_size_v<RowIndexes>,
//                                std::tuple_size_v<ColumnIndexes>>,
//                  RowIndexes, ColumnIndexes>;

// template <typename... Types>
// using column_vector =
//     typed_column_vector<Eigen::Vector<representation, sizeof...(Types)>,
//                         Types...>;

// template <typename... Types>
// using row_vector =
//     typed_row_vector<Eigen::RowVector<representation, sizeof...(Types)>,
//                      Types...>;

//! @brief Strongly typed linear algebra samples.
//!
//! @details A variety of activities of strongly typed linear algebra with
//! std::mdspan, std::linalg, and mp-units.
[[maybe_unused]] auto sample{[] { return 0; }()};
} // namespace
} // namespace sample
} // namespace fcarouge
